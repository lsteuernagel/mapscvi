
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
#' @param max_epochs epochs to train
#' @param assay assay name for new prediction and where to find raw counts. defaults to RNA
#' @param use_reticulate if TRUE: tries to call scvi via reticulate, if FALSE: exports an anndata object to a temp file and then calls a python script to run prediction code using system()
#' @param temp_dir directory for temporary files if use_reticulate == FALSE
#' @param global_seed seed
#'
#' @return query_seurat_object: Updated seurat object with predicted latent space as reduction.
#'
#' @export
#'
#' @import SeuratObject Seurat SeuratDisk
#'
#' @examples

# TODO: need to limit cores used by scvi ! (run setup !)

predict_query = function(query_seurat_object,model_path,query_reduction="scvi",max_epochs = 30,assay="RNA",use_reticulate = FALSE,temp_dir = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/",global_seed=12345){

  # load the variable feature from modelpath
  var_features = utils::read.table(paste0(model_path,"var_names.csv"),header = F)$V1

  # export to anndata
  var_df = data.frame(var_names = rownames(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay)))
  rownames(var_df) = var_df$var_names
  # make a matrix with included variable genes in query
  included_var_features = intersect(var_features,var_df$var_names)
  # if I don't subset here, then a bigger anndata will be exported but the library size prior will be estimated correctly (update: no because the original model only used the x var genes to estimate lib size)
  matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay)[included_var_features,])
  # for scvi to work we also need the others: we add them as all zero columns at the end
  missing_var_features = setdiff(var_features,var_df$var_names)
  add_matrix = matrix(data = 0,nrow = length(missing_var_features),ncol = ncol(matrix_for_anndata))
  rownames(add_matrix) = missing_var_features
  colnames(add_matrix) = colnames(matrix_for_anndata)
  matrix_for_anndata = rbind(matrix_for_anndata,add_matrix)
  # make new var_df
  var_df = data.frame(var_names = rownames(matrix_for_anndata))
  rownames(var_df) = var_df$var_names
  message("Matrix for anndata dim ",dim(matrix_for_anndata)[1]," ",dim(matrix_for_anndata)[2])

  ### reticulate code section:
  if(use_reticulate){

    if (!requireNamespace("reticulate", quietly = TRUE)) {
      warning("The reticulate package must be installed to use this function when use_reticulate is set to TRUE.")
      return(NULL)
    }

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
    # put raw data in X slot!
    # adata_query$X = adata_query$raw$X$copy()
    # load_query_data
    vae_q = scvi$model$SCVI$load_query_data(
      adata = adata_query,
      reference_model = model_path, # get model directly from disk!
      inplace_subset_query_vars=TRUE
    )
    # train
    message(max_epochs)
    vae_q$train(max_epochs=as.integer(max_epochs),
                plan_kwargs=list(weight_decay=0.0)
    ) # use same epochs and weight_decay = 0
    # results

    #adata_query$obsm["X_scVI"] = vae_q$get_latent_representation() # get laten dim
    # get results directly into R dataframe
    scvi_prediction = vae_q$get_latent_representation()
    scvi_prediction = as.matrix(scvi_prediction)
    colnames(scvi_prediction) = paste0("scVI_",1:ncol(scvi_prediction))
    rownames(scvi_prediction) = colnames(matrix_for_anndata)

  }else{ # this version does not use reticulate to execute scvi

    # make path if not existing
    system(paste0("mkdir -p ",temp_dir))

    # make Seurat from updated matrix
    temp_seurat = SeuratObject::CreateSeuratObject(counts = matrix_for_anndata, meta.data = query_seurat_object@meta.data, project = query_seurat_object@project.name)
    # export anndata
    h5Seurat_filename=paste0(temp_dir,"temp_",temp_seurat@project.name,".h5Seurat")
    SeuratDisk::SaveH5Seurat(object = temp_seurat,filename = h5Seurat_filename,overwrite = TRUE)
    updated_name = gsub(".h5Seurat",paste0("_",assay,".h5ad"),h5Seurat_filename)
    SeuratDisk::Convert(h5Seurat_filename, dest = updated_name,assay=assay,verbose=FALSE,overwrite=TRUE)

    # call python script
    output_file = paste0(temp_dir,"predicted_",temp_seurat@project.name,".txt")
    system(paste0("python3 -u ",system.file("python/map_scvi.py", package = "mapscvi")," ",updated_name," ",model_path," ",output_file," ",max_epochs))
    #system(paste0("python3 -u inst/python/map_scvi.py ",updated_name," ",model_path," ",output_file," ",max_epochs))
    # system(paste0("python3 -u python/map_scvi.py ",updated_name," ",model_path," ",output_file," ",max_epochs))
    #system.file("inst/python/map_scvi.py",package = "mapscvi",lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.0/mapscvi/")
    # load results into R
    scvi_prediction = utils::read.table(output_file)
    scvi_prediction = as.matrix(scvi_prediction)
    colnames(scvi_prediction) = paste0("scVI_",1:ncol(scvi_prediction))
    rownames(scvi_prediction) = colnames(matrix_for_anndata)

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
#' @param n_neighbors n neighbors passed to ProjectUMAP
#' @param annoy.metric 'cosine' or 'euclidean'
#' @param label_vec a vector with labels from reference that will be propagated to query (requires same order as query_reduction!). defaults to NULL (nothing will be propagated). See also 'propagate_labels'
#' @param global_seed seed
#'
#' @return query_seurat_object: Updated seurat object with projected UMAP, labels and knn graph
#'
#' @export
#'
#' @import Seurat SeuratObject
#'
#' @examples

project_query = function(query_seurat_object,reference_map_reduc,reference_map_umap,query_reduction="scvi",assay="RNA",n_neighbors = NULL,
                         annoy.metric = "cosine",label_vec =NULL,global_seed=12345){

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

  ### project query onto umap
  message(Sys.time(),": Project UMAP.." )

  if(is.null(n_neighbors)){n_neighbors=reference_map_umap@misc$model$n_neighbors}

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

  message(Sys.time(),": Add results to Seurat.." )
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

  # propagate labels
  if(!is.null(label_vec)){
    message(Sys.time(),": Predict labels..." )
    labels <- tryCatch({
      propagate_labels(nn_idx=query_seurat_object@neighbors[["query_ref_nn"]]@nn.idx,label_vec=label_vec)
    },
    error=function(cond) {
      message("Cannot propagate labels. Returning NA. Error:",cond)
      return(NA)
    })
  }
  # add to seurat
  query_seurat_object@meta.data[,paste0("predicted")] = labels
  #return
  return(query_seurat_object)
}

##########
### propagate_labels
##########

#' Propagate labels to query
#'
#' Use kNN to propagate labels (any metadata column from reference) to query.
#'
#' @param nn_idx kNN
#' @param label_vec vector with labels, ids must be consistent with NN indices in kNN
#'
#' @return predicted label for each cell (kNN row)
#'
#' @export
#'
#' @examples

propagate_labels = function(nn_idx,label_vec){
  # check
  if(max(nn_idx)>length(label_vec)){
    message("Please provide a vector with labels corresponding to nn_idx")
    stop()
  }
  # get Knn based labels
  knn_labels = apply(nn_idx,1,function(row,cell_labels){
    freq_label = table(cell_labels[row])
    names(freq_label)[freq_label == max(freq_label)][1]
  }, cell_labels = label_vec)

  return(knn_labels)
}

##########
### map_new_seurat
##########

#' Map a query seurat onto a reference
#'
#' Wrapper functions that executes low level functions to prepare, predict and project new data
#' NOT IMPLEMENTED YET!
#'
#' @inheritParams prepare_query
#' @inheritParams predict_query
#' @inheritParams project_query
#' @param label_col the column name in reference_map_metadata with labels to propagate
#' @param reference_seurat reference seurat
#' @param reference_reduction name of scvi reduction in reference. Also expect a umap named as 'umap_'reference_reduction
#'
#' @return formatted seurat object
#'
#' @export
#'
#' @import SeuratObject Seurat
#'
#' @examples

map_new_seurat_hypoMap = function(query_seurat_object,suffix="query",assay="RNA",subset_col="",label_col="",subset_values=NULL,max_epochs,reference_seurat=NULL,reference_reduction="scvi",model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/hypothalamus_neurons_reference_model/",
                                  use_reticulate = FALSE,temp_dir = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/",global_seed=12345){

  # prepare
  query_seurat_object = prepare_query_hypoMap(query_seurat_object,suffix=suffix,assay=assay,subset_col=subset_col,subset_values=subset_values,normalize=TRUE,
                                              covariates=c(batch_var = "Batch_ID",inferred_sex = "inferred_sex.x",rpl_signature_expr_median = "rpl_signature_expr_median"),global_seed=global_seed)
  # predict with scvi
  query_seurat_object = predict_query(query_seurat_object,model_path,max_epochs = max_epochs,assay="RNA",global_seed=global_seed)

  # check
  if(is.null(reference_seurat)){
    stop("Please provide reference seurat with latent space, umap and metadata")
  }
  if(! (reference_reduction %in% names(reference_seurat@reductions) & paste0("umap_",reference_reduction) %in% names(reference_seurat@reductions))){
    stop("Cannot find '",reference_reduction,"' or 'umap_",reference_reduction,"'  in provided reference_seurat.")
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
