
# This file contains functions to visualize projected data.

##########
### plot_query_labels
##########

#' Plot labels on query and reference to compare visually.
#'
#' Is a wrapper around Seurat::DimPlot functions that allows to plot qeury and reference side-by-side or on top of each other
#' The overlay parameters allow to change the appearance of query points and the alpha of the reference points.
#'
#' @param query_seura_object query seurat object
#' @param reference_seurat reference seurat object
#' @param label_col the column name in reference_map metadata with labels to propagate
#' @param label_col_query he column name in query_seura_object metadata with labels to propagate
#' @param overlay overlay query onto label
#' @param overlay_color color of overlay points
#' @param overlay_alpha = alpha for overlay plot
#' @param query_pt_size numeric for pt size of query cells. defaults to NULL which will inherit the pt.size from the reference DimPlot
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
#' @import SeuratObject Seurat cowplot ggplot2
#'
#' @examples
#'
#'

plot_query_labels = function(query_seura_object,reference_seurat,label_col,label_col_query = "predicted", overlay = FALSE, overlay_color = "red", overlay_alpha = 0.2,query_pt_size=NULL, query_umap = "umap_scvi",reference_umap="umap_scvi",labelonplot=TRUE,noaxes=TRUE,nolegend=TRUE,...){

  # check
  if(is.null(reference_seurat)){stop("Please provide reference seurat with latent space, umap and metadata")}
  if(! (reference_umap %in% names(reference_seurat@reductions))){stop("Cannot find '",reference_umap,"' in provided reference_seurat.") }
  if(! (query_umap %in% names(query_seura_object@reductions))){ stop("Cannot find '",reference_umap,"' in provided reference_seurat.")}
  if(! label_col %in% colnames(reference_seurat@meta.data)){stop("Cannot find '",label_col,"' in reference_seurat to label data.") }

  # overlay mode
  if(overlay){
    # extract data for overlay from query
    plot_data = cbind(query_seura_object@reductions[[query_umap]]@cell.embeddings,query_seura_object@meta.data)
    # plot reference UMAP
    p_full=DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap,label = labelonplot,...)
    if(noaxes){p_full = p_full+Seurat::NoAxes()}
    if(nolegend){ p_full = p_full+Seurat::NoLegend()}
    # adjust alpha
    p_full[[1]]$layers[[1]]$aes_params$alpha = min(overlay_alpha,1)
    # save and remove geom_text layer
    if(labelonplot){
      save_geom_text = p_full$layers[[2]]
      p_full$layers[[2]] =NULL
    }
    # get pt size
    if(is.null(query_pt_size)){
      pt_size = p_full[[1]]$layers[[1]]$aes_params$size
    }else{
      pt_size = query_pt_size
    }
    # plot query points on top
    p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2]),size=pt_size,color=overlay_color)
    # add geom_labels back
    if(labelonplot){p_full$layers[[3]] = save_geom_text}
  }else{
    # need labels in query
    if(! label_col_query %in% colnames(query_seura_object@meta.data)){
      stop("Cannot find '",label_col_query,"' in query_seura_object to label data.")
    }
    # browser()
    # side-by-side
    xlims = c(min(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,1])-0.5,max(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,1])+0.5)
    ylims = c(min(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,2])-0.5,max(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,2])+0.5)
    p1 = Seurat::DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap,label = labelonplot,...)+xlim(xlims)+ylim(ylims)
    # get pt size
    if(is.null(query_pt_size)){
      pt_size = p1[[1]]$layers[[1]]$aes_params$size
    }else{
      pt_size = query_pt_size
    }
    p2 = Seurat::DimPlot(query_seura_object,group.by = label_col_query,reduction = query_umap,label = labelonplot,pt.size = pt_size,...)+xlim(xlims)+ylim(ylims)
    if(noaxes){
      p1 = p1+Seurat::NoAxes()
      p2 = p2+Seurat::NoAxes()
    }
    if(nolegend){
      p1 = p1+Seurat::NoLegend()
      p2 = p2+Seurat::NoLegend()
    }
    p_full = cowplot::plot_grid(p1,p2)

  }

  p_full
}


##########
### plot propagation statistics
##########

#' Plot propagation statistics about query column
#'
#' Determines over or und represented clusters
#' TODO: add description
#'
#' @param query_seura_object query seurat object
#' @param reference_seurat reference seurat object
#' @param label_col the column name in reference_map metadata with labels to propagate
#' @param label_col_query he column name in query_seura_object metadata with labels to propagate
#' @param return_data return data and / or  plot
#'
#' @return plot or data behind plot
#'
#' @export
#'
#' @import SeuratObject Seurat cowplot ggplot2
#'
#' @examples
#'
#'

plot_propagation_stats = function(query_seura_object,reference_seurat,label_col,label_col_query = "predicted",return_data=FALSE){

  # TODO: add code!


  if(return_data){
    return(propagation_stats)
  }
}


##########
### plot_cluster_tree
##########

#' Plot metadata entries in cluster tree using ggtree
#'
#' TODO: add description
#' TODO: add reference to ggtree vignette
#'
#' @param treedata TODO
#' @param metadata TODO
#' @param leaf_level which level to use as leaves ?
#'
#' @return ggtree object or plot
#'
#' @export
#'
#' @import SeuratObject Seurat cowplot ggplot2
#'
#' @examples
#'
#'

plot_cluster_tree = function(treedata,metadata,leaf_level){

  # optional use of packages: ggtree imports also tidytree and treeio!
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    warning("The ggtree package must be installed to use this functionality")
    #Either exit or do something without rgl
    return(NULL)
  }

  # TODO: add code!

}



##########
### plot_sankey_comparison
##########

#' Plot comparison between reference and query using sankey plots
#'
#' TODO: add description
#'
#' @param metadata TODO
#' @param return_data
#'
#' @return ggtree object or plot
#'
#' @export
#'
#' @import SeuratObject Seurat cowplot ggplot2
#'
#' @examples
#'
#'

plot_sankey_comparison = function(metadata,return_data=FALSE){

  # optional use of packages:
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    warning("The networkD3 package must be installed to use this functionality")
    #Either exit or do something without rgl
    return(NULL)
  }
  #

  # TODO: add code!
  # make edges == edge_list_chord and nodes

  # p <- networkD3::sankeyNetwork(Links = edges, Nodes = nodes,
  #                               Source = "IDsource", Target = "IDtarget",
  #                               Value = "n", NodeID = "name",
  #                               sinksRight=FALSE,fontSize=20)
  # p

  if(return_data){
    return(list(nodes = nodes, edges = edges))
  }

}
