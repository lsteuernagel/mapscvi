
# This file contains the core functions to embed seurat data into a scvi space and an associated umap.
# These functions require execution of external python scripts!

##########
### predict_query
##########

#' Predict latent space of embedding
#'
#' Function to run scvi model prediction on query seurat based on provided scvi model.
#' Uses the scArches algorithm as described in the documentation: https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.load_query_data.html
#'
#' @param query_seurat_object Seurat object
#' @param model_path path to pretrained scvi model
#' @param query_reduction name for reduction
#' @param var_names optionally provide features to subset anndata to
#' @param max_epochs epochs to train
#' @param assay assay name for new prediction and where to find raw counts. defaults to RNA
#' @param use_reticulate if TRUE: tries to call scvi via reticulate, if FALSE: exports an anndata object to a temp file and then calls a python script to run prediction code using system()
#' @param pythonEnv path to a python binary, passed to reticulate::use_python
#' @param global_seed seed
#'
#' @return query_seurat_object: Updated seurat object with predicted latent space as reduction.
#'
#' @export
#'
#' @import SeuratObject Seurat SeuratDisk
#'
#' @importFrom data.table fread

# TODO: need to limit cores used by scvi ! (run setup !)

predict_query = function(query_seurat_object,model_path,query_reduction="scvi",var_names=NULL,max_epochs = 30,assay="RNA",pythonEnv = NULL,use_reticulate = FALSE,global_seed=12345){


  # check if model path is valid:
  if(!file.exists(paste0(model_path,"/","model.pt"))){
    stop("Error: Please provide a valid model_path to an scvi model at ",model_path)
  }

  # load the variable feature from modelpath
  if(!is.null(var_names)){
    # export to anndata
    var_df = data.frame(var_names = rownames(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay)))
    rownames(var_df) = var_df$var_names
    # make a matrix with included variable genes in query
    included_var_features = intersect(var_features,var_df$var_names)
    # if I don't subset here, then a bigger anndata will be exported but the library size prior will be estimated correctly (update: no because the original model only used the x var genes to estimate lib size)
    matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay)[included_var_features,])
    # make new var_df
    var_df = data.frame(var_names = rownames(matrix_for_anndata))
    rownames(var_df) = var_df$var_names
    message("Matrix for anndata dim ",dim(matrix_for_anndata)[1]," ",dim(matrix_for_anndata)[2])
  }else{
    matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay))
    # NEW CHANGE: make var_df object neccessary for sc$AnnData setup
    var_df = data.frame(var_names = rownames(matrix_for_anndata))
    rownames(var_df) = var_df$var_names
    message("Matrix for anndata dim ",dim(matrix_for_anndata)[1]," ",dim(matrix_for_anndata)[2])
  }

  ### reticulate code section:
  if(use_reticulate){

    if (!requireNamespace("reticulate", quietly = TRUE)) {
      warning("The reticulate package must be installed to use this function when use_reticulate is set to TRUE.")
      return(NULL)
    }

    # specifiy pythonEnv
    if (!is.null(pythonEnv)){
      reticulate::use_python(pythonEnv)
    }

    # TODO: add checks for python modules like this:
    #  py_module_available("scipy")
    # see also: https://rstudio.github.io/reticulate/articles/package.html

    # I am following this guide: https://docs.scvi-tools.org/en/stable/user_guide/notebooks/scvi_in_R.html
    pd <- reticulate::import('pandas', convert = FALSE)
    sc <- reticulate::import('scanpy', convert = FALSE)
    scvi <- reticulate::import('scvi', convert = FALSE)
    #scvi$settings$progress_bar_style = 'tqdm'

    # make anndata in python
    adata_query <- sc$AnnData(
      X   = t(matrix_for_anndata), #scVI requires raw counts, scanpy transposed data
      obs = query_seurat_object@meta.data,
      var = var_df
    )
    # put raw data in X slot ??

    # prepare:
    scvi$model$SCVI$prepare_query_anndata(adata_query, model_path)

    # adata_query$X = adata_query$raw$X$copy()
    # load_query_data
    vae_q = scvi$model$SCVI$load_query_data(
      adata = adata_query,
      reference_model = model_path
    )
    # train
    message(max_epochs)
    vae_q$train(max_epochs=as.integer(max_epochs),
                plan_kwargs=list(weight_decay=0.0),
                progress_bar_refresh_rate=0
    ) # use same epochs and weight_decay = 0

    # results
    # get results directly into R dataframe
    scvi_prediction = vae_q$get_latent_representation()
    scvi_prediction = as.matrix(scvi_prediction)
    colnames(scvi_prediction) = paste0("scVI_",1:ncol(scvi_prediction))
    rownames(scvi_prediction) = colnames(matrix_for_anndata)

  }else{ # this version does not use reticulate to execute scvi

    # tempfiles - should be the same across the session
    temp_dir = paste0(tempdir(),"/")

    # make Seurat from updated matrix
    temp_seurat = SeuratObject::CreateSeuratObject(counts = matrix_for_anndata, assay = assay, meta.data = query_seurat_object@meta.data, project = query_seurat_object@project.name)
    # export anndata
    h5Seurat_filename=paste0(temp_dir,"temp_",temp_seurat@project.name,".h5Seurat")
    temp_seurat[[assay]] <- as(object = temp_seurat[[assay]], Class = "Assay") # change assay from assay5 to fix 'Error in guess_dtype(x = x, ...) : unknown type'
    SeuratDisk::SaveH5Seurat(object = temp_seurat,filename = h5Seurat_filename,overwrite = TRUE)
    updated_name = gsub(".h5Seurat",paste0("_",assay,".h5ad"),h5Seurat_filename)
    SeuratDisk::Convert(h5Seurat_filename, dest = updated_name,assay=assay,verbose=FALSE,overwrite=TRUE)

    # call python script
    output_file = paste0(temp_dir,"predicted_",temp_seurat@project.name,".txt")
    system(paste0("python3 -u ",system.file("python/map_scvi2.py", package = "mapscvi")," ",updated_name," ",model_path," ",output_file," ",max_epochs))
    #system(paste0("python3 -u inst/python/map_scvi.py ",updated_name," ",model_path," ",output_file," ",max_epochs))
    # system(paste0("python3 -u python/map_scvi.py ",updated_name," ",model_path," ",output_file," ",max_epochs))
    #system.file("inst/python/map_scvi.py",package = "mapscvi",lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.0/mapscvi/")
    # load results into R
    scvi_prediction = data.table::fread(output_file,header = TRUE,data.table = F)
    #scvi_prediction = as.matrix(scvi_prediction)
    rownames_x = as.character(scvi_prediction[,1])
    #print(paste0(rownames_x,collapse = " | "))
    scvi_prediction = scvi_prediction[,2:ncol(scvi_prediction)]
    # sensure numeric
    #scvi_prediction[,apply(scvi_prediction,2,is.character)] <- apply(scvi_prediction[,apply(scvi_prediction,2,is.character)],2, as.numeric)
    scvi_prediction=as.matrix(apply(scvi_prediction,2,as.numeric))
    rownames(scvi_prediction) = rownames_x
    #print(paste0(apply(scvi_prediction,2,is.numeric),collapse = " | "))
    # colnames(scvi_prediction) = paste0("scVI_",1:ncol(scvi_prediction))


    # scvi_prediction = ...
  }

  # make query dim red
  query_dimred <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(scvi_prediction),
    stdev = as.numeric(apply(scvi_prediction, 2, stats::sd)),
    assay = assay,
    key = query_reduction
  )
  # add to query seurat
  query_seurat_object@reductions[[query_reduction]] = query_dimred

  # return result
  return(query_seurat_object)
}

##########
### project_query
##########

#' Project query into reference UMAP
#'
#' Function to project query onto reference using a predicted scvi latent space. Returns UMAP, propagated labels and kNN in the updated seurat object.
#'
#' @param query_seurat_object Seurat object
#' @param reference_map_reduc reference map dim_reduc (scvi latent space)
#' @param reference_map_umap reference map UMAP dim_reduc
#' @param query_reduction key of predicted reduction used for UMAP projection
#' @param assay assay name. defaults to RNA
#' @param use_projectUMAP whether to use Seurat's projectUMAP  or a manual adaption (circumventing possible problems between Seurat and uwot)
#' @param n_neighbors n neighbors passed to ProjectUMAP
#' @param annoy.metric 'cosine' or 'euclidean'
#' @param result_type return type: "all","label","entropy" .  For label the maximum is taken to return the likeliest assignment. All returns all results (multiple colummns!), entropy calculates the entropy across probabilities as an uncertainty measure but does not return a label
#' @param label_vec a vector with labels from reference that will be propagated to query (requires same order as query_reduction!). defaults to NULL (nothing will be propagated). See also 'propagate_labels'
#' @param global_seed seed
#'
#' @return query_seurat_object: Updated seurat object with projected UMAP, labels and knn graph
#'
#' @export
#'
#' @import Seurat SeuratObject uwot

project_query = function(query_seurat_object,reference_map_reduc,reference_map_umap,query_reduction="scvi",assay="RNA",use_projectUMAP=FALSE,
                         annoy.metric = "cosine",result_type="label",label_vec =NULL,global_seed=12345){

  # TODO: test for umap model (?)

  #test for query_reduction name
  if(! query_reduction %in% names(query_seurat_object@reductions)){
    message("Error: query_reduction name cannot be found in query_seurat_object")
    stop()
  }
  #test for same space size
  if(ncol(query_seurat_object@reductions[[query_reduction]]@cell.embeddings) != ncol(reference_map_reduc@cell.embeddings)){
    message("Warning: query reduction and reference reduction don't have the same number of dimensions!")
    # stop()
  }

  #if(is.null(n_neighbors)){n_neighbors=reference_map_umap@misc$model$n_neighbors}
  n_neighbors=reference_map_umap@misc$model$n_neighbors

  if(use_projectUMAP){

    ### project query onto umap
    message(Sys.time(),": Project UMAP using Seurat::ProjectUMAP.." )

    query_projection = Seurat::ProjectUMAP(
      query = query_seurat_object@reductions[[query_reduction]],
      query.dims = 1:ncol(query_seurat_object@reductions[[query_reduction]]@cell.embeddings),
      reference = reference_map_reduc,
      reference.dims = 1:ncol(reference_map_reduc@cell.embeddings),
      k.param = reference_map_umap@misc$model$n_neighbors,
      n.neighbors = n_neighbors,
      nn.method = "annoy",
      n.trees = 50,
      annoy.metric = annoy.metric,
      l2.norm = FALSE,
      seed_use = global_seed,
      neighbor.name = "query_ref_nn",
      reduction.model = reference_map_umap,
      return.neighbor = TRUE
    )
    # make a new dim red (just to be sure that key and colnames are coherent)
    query_umap <- Seurat::CreateDimReducObject(
      embeddings = query_projection$proj.umap@cell.embeddings,
      loadings = query_projection$proj.umap@feature.loadings,
      projected = query_projection$proj.umap@feature.loadings.projected,
      stdev = query_projection$proj.umap@stdev,
      assay = assay,
      key = paste0("umap_",query_reduction),
      jackstraw = query_projection$proj.umap@jackstraw
    )
    # add to query seurat
    query_seurat_object@reductions[[paste0("umap_",query_reduction)]] = query_umap

    # also add nn object to seurat
    query_seurat_object@neighbors[["query_ref_nn"]]=query_projection$query.neighbor

  }else{

    message(Sys.time(),": Project UMAP using manual uwot_transform.." )

    # extract model from umap reduc
    model <- Misc(
      object = reference_map_umap,
      slot = "model"
    )
    # we need to set num_precomputed_nns = 1 to avoid: Error in umap_transform(X = NULL, nn_method = query.neighbor, model = model,  :  Expecting
    model$num_precomputed_nns = 1

    message(Sys.time(),": Find neighbors in reference.." )

    # run manual neighbor mapping
    query.neighbor <- FindNeighbors(
      object = reference_map_reduc@cell.embeddings,
      query = query_seurat_object@reductions[[query_reduction]]@cell.embeddings,
      k.param = n_neighbors,
      nn.method = "annoy",
      n.trees = 50,
      annoy.metric = annoy.metric,
      cache.index = FALSE,
      index = NULL,
      return.neighbor = TRUE,
      l2.norm = FALSE
    )

    # make a uwot conformable object
    query.neighbor.uwot = list(idx = query.neighbor@nn.idx, dist = query.neighbor@nn.dist)

    message(Sys.time(),": Run umap_transform.." )

    # this way we can provide a precomputed nn object:
    transformed_query = uwot::umap_transform(
      X = NULL,
      nn_method = query.neighbor.uwot,
      model = model,
      n_threads = 30,
      n_epochs = NULL,
      verbose = F
    )

    message(Sys.time(),": Add results to Seurat.." )

    # make query dim red
    query_umap <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(transformed_query),
      stdev = as.numeric(apply(transformed_query, 2, stats::sd)),
      assay = assay,
      key = paste0("umap_",query_reduction),
    )

    # update cell ids:
    rownames(query_umap@cell.embeddings) = rownames(query_seurat_object@meta.data)

    # add to query seurat
    query_seurat_object@reductions[[paste0("umap_",query_reduction)]] = query_umap

    # also add nn object to seurat
    query_seurat_object@neighbors[["query_ref_nn"]]= query.neighbor

  }

  ##

  # propagate labels
  if(!is.null(label_vec)){
    if(nrow(reference_map_reduc@cell.embeddings) != length(label_vec)){
      message("Warning: label vector does not have the same length as reference reduction. Skipping prediction.")
      # stop()
    }else{
      message(Sys.time(),": Predict labels..." )
      # run label propagation
      query_seurat_object =propagate_labels_prob(neighbors_object=query_seurat_object@neighbors[["query_ref_nn"]],
                                                 query_seurat_object = query_seurat_object,
                                                 label_vec = label_vec,
                                                 with_euclidean =FALSE,
                                                 result_type = result_type,
                                                 add_to_seurat =TRUE)
      # additionally run neighbor distances qc
      query_seurat_object = avg_neighbor_distances(query_seurat_object = query_seurat_object,
                                                   reference_map_reduc = reference_map_reduc, # reference_map_reduc
                                                   query_nn = "query_ref_nn",
                                                   distance.metric=annoy.metric,
                                                   add_to_seurat=TRUE)
    }
  }

  #return
  return(query_seurat_object)
}

##########
### adjusted_cell_probabilities
##########

#' Calculate probabilities for labels of neighboring cells
#'
#' Cell probabilities similar to scARches algorithm
#' This function runs per cell.
#'
#' @param dist_Nc vector of length k with euclidean distances to neigbors
#' @param labels_of_neighbor vector of length k with labels of neighbors
#' @param with_euclidean if input is euclidean: use gussian smoothing and convert to similarity through this. If not assume cosine distance and just invert using: 1 - dist_Nc^2/2
#' @param result_type return type: "all","label","entropy" .  For label the maximum is taken to return the likeliest assignment. All returns all results (multiple colummns!), entropy calculates the entropy across probabilities as an uncertainty measure but does not return a label
#'
#' @return label probability depending on result_type
#'
#' @export

adjusted_cell_probabilities = function(dist_Nc,labels_of_neighbor,with_euclidean = FALSE,result_type="all"){

  # step 1: Distances and input
  # ... input of function
  local_labels = unique(labels_of_neighbor)
  k = length(labels_of_neighbor)
  # apply gaussian ?
  if(with_euclidean){
    # step 2: compute the standard deviation of the nearest distances
    sd_nc = sqrt(sum(dist_Nc^2) / k )
    # step 3: apply Gaussian kernel to distances
    d_app = exp(-1*(dist_Nc/(2/sd_nc)^2) )
  }else{
    d_app = 1 - dist_Nc^2/2 # else just invert cosine distance
  }
  # step 4: we computed the probability of assigning each label y to the query cell c by normalizing across all adjusted distances
  label_probabilities = tapply(d_app,INDEX = labels_of_neighbor,FUN = sum) / sum(d_app)
  # return result
  if(result_type == "all"){
    return(label_probabilities)
    # or summarise further before returning:
  }else{
    if(result_type == "label"){
      # step 5: uncertainty score: 1 - prob_label --> I rather return directly the highest label with its name
      return(label_probabilities[label_probabilities == max(label_probabilities)][1])
    }else if(result_type == "entropy"){
      # step 6: try out entropy:
      ## entropy helper:
      entropy_fun = function(x,logfun ="log2"){
        log_vec = do.call(logfun,list(x))
        log_vec[is.infinite(log_vec)] = 0
        log_vec[is.nan(log_vec)] = 0
        return(-sum(x * log_vec))
      }
      entropy_uncertainty = entropy_fun(label_probabilities) / max(1,log2(length(local_labels)))
      return(entropy_uncertainty)
    }else{
      return(NA)
    }
  }

}

##########
### propagate_labels_prob
##########

#' Propagates cell labels base don provided vector using the nearest neighbors in the reference.
#'
#' Cell probabilities similar to scARches algorithm
#' This function runs per cell
#'
#' @param neighbors_object neighbors-object with reference neighbors of query (see SeuratObject). If NULL will run neighbor detection between query_seurat_object and reference_seurat_object. Requires slot nn.idx and nn.dist.
#' @param label_vec vector with labels
#' @param query_seurat_object seurat with query data (and a reduction that is projected from reference_seurat_object) to find shared neighbors. NULL by default (not used when neighbors_object is provided)
#' @param reference_seurat_object seurat to map onto. NULL by default (not used when neighbors_object is provided)
#' @param reduction_name_query name of reduction in query (not used when neighbors_object is provided)
#' @param reduction_name_reference name of reduction in reference (not used when neighbors_object is provided)
#' @param annoy.metric euclidean or cosine (not used when neighbors_object is provided)
#' @param k.param k param for neighbor finding (not used when neighbors_object is provided)
#' @param with_euclidean Apply gaussian kernel to smooth distances .only use when distance = euclidean ! (or a neighbors_object based on euclidean distances is provided)
#' @param result_type "all","label" or "entropy". All means all probabilities of cell types in neighborhood. For label the maximum is taken to return the likeliest assignment. All returns all results (multiple colummns!), entropy calculates the entropy across probabilities as an uncertainty measure but does not return a label
#' @param add_entropy re-run to calculate entropy as uncertainty measure
#' @param add_to_seurat add to query_seurat_object or return dataframe, requires query_seurat_object to be provided
#'
#' @return seuratobject or dataframe with label propagation results and qc
#'
#' @importFrom dplyr bind_rows
#'
#' @export

propagate_labels_prob = function(neighbors_object=NULL,label_vec,query_seurat_object=NULL,reference_seurat_object=NULL,reduction_name_query="scvi",reduction_name_reference="scvi",annoy.metric="cosine",k.param=30, with_euclidean =FALSE,result_type = "label", add_entropy =FALSE, add_to_seurat =FALSE){

  if(is.null(label_vec)){stop("Error: Please provide a valid (non NULL) label_vec.")}

  # need euclidean distances neighbors
  if(is.null(neighbors_object)){
    if(ncol(reference_seurat_object) != length(label_vec)){
      stop("Error: Please provide a label_vec of the same length a cells in reference_seurat_object when calculating a new neighbor object.")
      # stop()
    }
    message("Calculating neighbor object")
    neighbors_object = Seurat::FindNeighbors(reference_seurat_object@reductions[[reduction_name_reference]]@cell.embeddings,
                                             query = query_seurat_object@reductions[[reduction_name_query]]@cell.embeddings,
                                             k.param = k.param, return.neighbor =TRUE,
                                             annoy.metric=annoy.metric)
  }else{
    if(class(neighbors_object)[1] != "Neighbor"){
      stop("Error: Please provide a valid neighbors_object.")
    }
    if(max(neighbors_object@nn.idx) > length(label_vec)){
      stop("Error: Please provide a label_vec of the same length as reference indices in neighbors_object@nn.idx.")
    }
  }
  # apply max prob per cell function
  if(with_euclidean & annoy.metric=="cosine"){message("Warning: Applying gaussian filter after using cosine distance.")}
  message("Estimate probabilities")
  n=nrow(neighbors_object@nn.dist)
  probabilities = sapply(1:n,function(x,distances,neighbor_idxs,labels){
    dist_Nc = distances[x,]
    label_of_neighbor = labels[neighbor_idxs[x,]]
    prob = adjusted_cell_probabilities(dist_Nc = dist_Nc,labels_of_neighbor = label_of_neighbor,with_euclidean = with_euclidean,result_type = result_type)
    prob
  },distances = neighbors_object@nn.dist,neighbor_idxs = neighbors_object@nn.idx,labels = label_vec)

  # mapping results
  mapping_results = data.frame(Cell_ID = neighbors_object@cell.names[1:n])
  # rename
  if(result_type == "label"){
    mapping_results = as.data.frame(cbind(mapping_results,as.numeric(probabilities)))
    colnames(mapping_results)[ncol(mapping_results)] = "prediction_probability"
    mapping_results$predicted = names(probabilities)
  }
  if(result_type == "entropy"){
    mapping_results = as.data.frame(cbind(mapping_results,as.numeric(probabilities)))
    colnames(mapping_results)[ncol(mapping_results)] = "prediction_entropy"
  }
  if(result_type == "all"){
    probabilities_df <- data.frame(do.call(bind_rows, probabilities))
    probabilities_df[is.na(probabilities_df)] = 0
    mapping_results = as.data.frame(cbind(mapping_results,probabilities_df))
    colnames(mapping_results)[2:ncol(mapping_results)] = paste0("prediction_",colnames(probabilities_df))
  }
  # return
  if(add_to_seurat & !is.null(query_seurat_object)){
    # remove if existing for clean join
    keep_names = colnames(query_seurat_object@meta.data)[!colnames(query_seurat_object@meta.data) %in% c("predicted","prediction_probability","prediction_entropy")]
    keep_names = keep_names[ ! grepl("prediction_",keep_names)]
    query_seurat_object@meta.data = query_seurat_object@meta.data[,keep_names]
    # join
    query_seurat_object@meta.data = dplyr::left_join(query_seurat_object@meta.data, mapping_results, by = c("Cell_ID"="Cell_ID"))
    # set rownames of dataframe after dplyr
    rownames(query_seurat_object@meta.data) = query_seurat_object@meta.data$Cell_ID
    return(query_seurat_object)
  }else{
    return(mapping_results)
  }
}


##########
### map_new_seurat_hypoMap
##########

#' Map a query seurat onto the hypomap reference
#'
#' Wrapper functions that executes low level functions to prepare, predict and project new data onto a reference based on an scvi model.
#' Defaults to a smaller object (that can be shipped with the package) for reference_seurat (mapscvi::reference_hypoMap_downsample). Please alternatively download and specify the full object: https://www.repository.cam.ac.uk/handle/1810/341048
#' If the model differs strongly from the hypoMap model (covariates used etc) the prepare and mapping functions might break.
#' In that case please manually run the prepare, predict and project functions as shown in the readme.
#'
#' @inheritParams prepare_query
#' @inheritParams predict_query
#' @inheritParams project_query
#' @param label_col the column name in reference_map_metadata with labels to propagate
#' @param reference_mode directly specify what to use as reference. Possible values are either 'hypoMap_neurons', 'hypoMap_full' or 'manual'. This will set default values for reference_seurat, reference_reduction and model_path that can be overwritten by directly specifying them. Setting this parameter to 'manual' or NULL requires valid entries for the other parameters.
#' @param reference_seurat reference seurat
#' @param reference_reduction name of scvi reduction in reference. Also expect a umap named as 'umap_'reference_reduction
#' @param inferred_sex_varname variable name for inferred sex that has to correspond to model (will automatically be set when using reference_mode)
#'
#' @return formatted seurat object
#'
#' @export
#'
#' @import SeuratObject Seurat

map_new_seurat_hypoMap = function(query_seurat_object,suffix="query",assay="RNA",subset_col="",
                                  label_col="",subset_values=NULL,max_epochs,
                                  model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", package = 'mapscvi'),
                                  reference_seurat=mapscvi::reference_hypoMap_downsample,reference_reduction="scvi",inferred_sex_varname = "inferred_sex" ,
                                  pythonEnv = NULL, use_reticulate = FALSE,global_seed=12345){

  # prepare
  query_seurat_object = prepare_query(query_seurat_object,suffix=suffix,assay=assay,subset_col=subset_col,subset_values=subset_values,normalize=TRUE,batch_var = "Batch_ID",global_seed=global_seed)

  # check if model path is valid:
  if(!file.exists(paste0(model_path,"/","model.pt"))){
    stop("Error: Please provide a valid model_path to an scvi model at ",model_path)
  }
  # predict with scvi
  query_seurat_object = predict_query(query_seurat_object,
                                      model_path,
                                      max_epochs = max_epochs,
                                      assay="RNA",
                                      use_reticulate = use_reticulate,
                                      pythonEnv = pythonEnv,
                                      global_seed=global_seed)
  # check
  if(is.null(reference_seurat)){
    stop("Error: Please provide reference seurat with latent space, umap and metadata")
  }
  if(! (reference_reduction %in% names(reference_seurat@reductions) & paste0("umap_",reference_reduction) %in% names(reference_seurat@reductions))){
    stop("Error: Cannot find '",reference_reduction,"' or 'umap_",reference_reduction,"'  in provided reference_seurat.")
  }

  # make label_vec
  reference_map_metadata = reference_seurat@meta.data
  label_vec=NULL
  if(!is.null(reference_map_metadata)){
    if(label_col %in% colnames(reference_map_metadata)){
      label_vec = reference_map_metadata[,label_col]
    }else{
      message("label_col does not exist in provided metadata. Cannot propagate labels!")
    }
  }else{
    message("No metadata provided. Cannot propagate labels!")
  }
  # project onto reference
  query_seurat_object = project_query(query_seurat_object,reference_seurat@reductions[[reference_reduction]],
                                      reference_seurat@reductions[[paste0("umap_",reference_reduction)]],
                                      label_vec =label_vec,
                                      global_seed=global_seed)

  # return
  return(query_seurat_object)
}
