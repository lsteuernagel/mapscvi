
# This file contains functions to visualize projected data.

plot_query_labels = function(query_seura_object,label_col,label_col_query = NULL, overlay = FALSE, query_reduction = paste0("umap_","scvi"),reference_seurat=NULL,reference_map_umap=NULL,reference_map_metadata=NULL,dim_red_file){

  # gather data
  if(!is.null(reference_seurat)){
    if(query_reduction %in% names(reference_seurat@reductions)){
      reference_map_umap = reference_seurat@reductions[[query_reduction]]
    }else{
      if(is.null(reference_map_umap)){
        stop("Cannot find ",query_reduction," in Seurat object. Stop.")
      }else{
       # message("Cannot find ",query_reduction," in Seurat object. Using provided reference_map_umap for plotting")
      }
    }
    reference_map_metadata = reference_seurat@meta.data
  }
  if(is.null(reference_map_umap) | is.null(reference_map_metadata)){
    stop("Please provide reference umap and metadata")
    # message("Loading reference information from ",dim_red_file)
    # load(dim_red_file)
  }

  # side-by-side

  # overlay mode

}
