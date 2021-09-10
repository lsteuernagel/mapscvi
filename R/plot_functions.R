
# This file contains functions to visualize projected data.

##########
### plot_query_labels
##########

#' Plot labels on query and reference to compare visually.
#'
#' Is a wrapper around DimPlot functions
#'
#' @inheritParams project_query
#' @param label_col the column name in reference_map metadata with labels to propagate
#' @param label_col_query he column name in query_seura_object metadata with labels to propagate
#' @param overlay overlay query onto label
#' @param overlay_color color of overlay points
#' @param overlay_alpha = alpha for overlay plot
#' @param query_umap name of query umap. defaults to "umap_scvi"
#' @param reference_umap name of reference umap. defaults to "umap_scvi"
#' @param labelonplot put labels on plot. defaults to TRUE
#' @param noaxes don't plot axes on UMAP. defaults to TRUE
#' @param nolegend don't plot legend. defaults to TRUE
#' @param ... additional arguments to Seurat::DimPlot
#'
#' @return plot
#'
#' @export
#'
#' @import SeuratObject Seurat
#'
#' @examples
#'
#'

plot_query_labels = function(query_seura_object,reference_seurat,label_col,label_col_query = NULL, overlay = FALSE, overlay_color = "orange", overlay_alpha = 0.2, query_umap = "umap_scvi",reference_umap="umap_scvi",labelonplot=TRUE,noaxes=TRUE,nolegend=TRUE,...){

  # check
  if(is.null(reference_seurat)){
    stop("Please provide reference seurat with latent space, umap and metadata")
  }
  if(! (reference_reduction %in% names(reference_seurat@reductions) & paste0("umap_",reference_reduction) %in% names(reference_seurat@reductions))){
    stop("Cannot find '",reference_reduction,"' or 'umap_",reference_reduction,"'  in provided reference_seurat.")
  }
  if(is.null(label_col_query)){
    if(label_col %in% colnames(query_seura_object@meta.data)){
      label_col_query = paste0(label_col
    }
    stop("Cannot find predicted_'",label_col,"' query_seura_object to label data. Please try to provide directly via the label_col_query argument.")
  }

  # overlay mode
  if(overlay){
    # extract data for overlay from query
    plot_data = cbind(query_seura_object@reductions[[query_umap]]@cell.embeddings,query_seura_object@meta.data)
    # plot reference UMAP
    p_full=DimPlot(reference_seurat,group.by = label_col,label.size = 5,reduction = reference_umap,label = labelonplot,...)
    if(noaxes){p_full = p_full+Seurat::NoAxes()}
    if(nolegend){ p_full = p_full+Seurat::NoLegend()}
    # adjust alpha
    p_full[[1]]$layers[[1]]$aes_params$alpha = min(overlay_alpha,1)
    # save and remove geom_text layer
    save_geom_text = p_full$layers[[2]]
    p_full$layers[[2]] =NULL
    # get pt size
    pt_size = p_full[[1]]$layers[[1]]$aes_params$size
    # plot query points on top
    p_full=p_full+geom_point(data=plot_data,aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2]),size=pt_size,color=overlay_color)
    # add geom_labels back
    if(labelonplot){p_full$layers[[3]] = save_geom_text}
  }else{
    # side-by-side
    p1 = Seurat::DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap,label = labelonplot,...)
    p1 = Seurat::DimPlot(query_seura_object,group.by = label_col_query,reduction = query_umap,label = labelonplot,...)
    if(noaxes){
      p1 = p1+Seurat::NoAxes()
      p2 = p2+Seurat::NoAxes()
    }
    if(nolegend){
      p1 = p1+Seurat::NoLegend()
      p2 = p2+Seurat::NoLegend()
    }
    p_full = cowplot::plot_grid(p1,p2)

    p_full=DimPlot(neuron_map_seurat,group.by = "K169_named",label.size = 1.5,label=TRUE)+Seurat::NoLegend()

  }

  p_full
}
