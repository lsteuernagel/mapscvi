
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
#' @param global_seed seed
#' 
#' @return query_seurat_object: Updated seurat object with predicted latent space as reduction.
#' 
#' @export
#' 
#' @import Seurat
#' @import SeuratObject
#' 
#' @examples

# TODO: need to limit cores used by scvi ! (run setup !)

predict_query = function(query_seurat_object,model_path,query_reduction="scvi",max_epochs = 30,assay="RNA",global_seed=12345){
  
  # load the variable feature from modelpath
  var_features = utils::read.table(paste0(model_path,"var_names.csv"),header = F)$V1
  
  # I am following this guide: https://docs.scvi-tools.org/en/stable/user_guide/notebooks/scvi_in_R.html 
  pd <- import('pandas', convert = FALSE)
  sc <- import('scanpy', convert = FALSE)
  scvi <- import('scvi', convert = FALSE)
  #scvi$settings$progress_bar_style = 'tqdm'
  
  # export to anndata
  var_df = data.frame(var_names = rownames(SeuratObject::GetAssayData(query_seurat_object,slot='counts',assay=assay)))
  rownames(var_df) = var_df$var_names
  # make a matrix with included variable genes in query
  included_var_features = intersect(var_features,var_df$var_names)
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
  # make dataframe with results (python version, need to fix some reticulate stuff):
  # output = pd$DataFrame(adata_query$obsm["X_scVI"])
  # output = output$set_index(adata_query$obs_names)
  # scvi_prediction = output$set_axis(("scVI_" + str(s) for s in output$axes[1]$to_list()), axis=1, inplace=False)
  # scvi_prediction = as.data.frame(scvi_prediction) # just to make sure we are back in R
  
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
#' @param annoy.metric 'cosine' or 'euclidean'
#' @param label_vec a vector with labels from reference that will be propagated to query (requires same order as query_reduction!). defaults to NULL (nothing will be propagated). See also 'propagate_labels'
#' @param global_seed seed
#' 
#' @return query_seurat_object: Updated seurat object with projected UMAP, labels and knn graph
#' 
#' @export
#' 
#' @import Seurat
#' @import SeuratObject
#' 
#' @examples

project_query = function(query_seurat_object,reference_map_reduc,reference_map_umap,query_reduction="scvi",assay="RNA",k_param_umap = 30,
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
  
  query_projection = Seurat::ProjectUMAP(
    query = query_seurat_object@reductions[[query_reduction]],
    query.dims = 1:ncol(query_seurat_object@reductions[[query_reduction]]@cell.embeddings),
    reference = reference_map_reduc,
    reference.dims = 1:ncol(reference_map_reduc@cell.embeddings),
    k.param = k_param_umap,
    n.neighbors = k_param_umap,
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
  if(!is.null(label_vec))
    message(Sys.time(),": Predict labels..." )
  labels <- tryCatch({
    propagate_labels(nn_idx=query_seurat_object@neighbors[["query_ref_nn"]]@nn.idx,label_vec=label_vec)
  },
  error=function(cond) {
    message("Cannot propagate labels. Returning NA. Error:",cond)
    return(NA)
  })
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
#' @param kNN kNN
#' @param label_col vector with labels, ids must be consistent with NN indices in kNN
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
#' 
#' @param query_seurat_object Seurat object
#' @param suffix project name to clearly label various steps. defaults to 'query'
#' @param subset_col column to subset by
#' @param subset_values values of subset_col to filter on
#' @param max_epochs epochs to train during scvi query function
#' @param model_path path to scvi model on disk
#' @param global_seed seed
#' 
#' @return formatted seurat object
#' 
#' @export
#' 
#' @examples

map_new_seurat_hypoMap = function(query_seurat_object,suffix="query",subset_col="",label_col="",subset_values=NULL,max_epochs,reference_map_reduc=NULL,reference_map_umap=NULL,reference_map_metadata=NULL,model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/hypothalamus_neurons_reference_model/",global_seed=12345){
  
  # prepare
  query_seurat_object = prepare_query_hypoMap(query_seurat_object,suffix="query",assay="RNA",subset_col=subset_col,subset_values=subset_values,normalize=TRUE,
                                              covariates=c(batch_var = "Batch_ID",inferred_sex = "inferred_sex.x",rpl_signature_expr_median = "rpl_signature_expr_median"),global_seed=global_seed)
  # predict with scvi
  query_seurat_object = predict_query(query_seurat_object,model_path,max_epochs = max_epochs,assay="RNA",global_seed=global_seed)
  
  # load reference map data
  if(is.null(reference_map_reduc) | is.null(reference_map_umap)){
    # TODO: load from rds!
    # reference_map_metadata!
  }
  # make label_vec
  label_vec=NULL
  if(!is.null(reference_map_metadata)){
    if(label_col %in% colnames(reference_map_metadata)){
      label_vec = reference_map_metadata[,label_col]
    }
  }
  # project onto reference
  query_seurat_object = project_query(query_seurat_object,reference_map_reduc,reference_map_umap,k_param_umap = 30,label_vec =label_vec,global_seed=global_seed)
  
  # return
  return(query_seurat_object)
}

### Test code:
# global_seed = 123467# seed
# map_name = "hypothalamus_neurons_reference"
# map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/"
# map_seurat_path = paste0(map_path,map_name,".h5Seurat")
# neuron_map_seurat = SeuratDisk::LoadH5Seurat(map_seurat_path)
# reference_map_reduc = neuron_map_seurat@reductions$scvi
# reference_map_umap = neuron_map_seurat@reductions$umap_scvi
# label_vec  = neuron_map_seurat@meta.data$K169_named
# 
# query_seurat_object_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/mapped_data_yeo_full/mapped_data_yeo_full.h5Seurat" # seurat object to load
# query_seurat_object = SeuratDisk::LoadH5Seurat(query_seurat_object_path)
