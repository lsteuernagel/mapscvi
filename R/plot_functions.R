
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
#' @param label_col_query he column name in query_seura_object metadata with labels to propagate
#' @param label_col the column name in reference_map metadata with labels to propagate
#' @param ratio_cut what ratio is considered to be a relevant difference (just for coloring)
#' @param text_size text_size in ggplots
#' @param plot_abs use absolute cell counts in barplot (might be option for large query data). defaults to FALSE
#' @param n_top how many clusters to include in barplots
#' @param ref_col color in plot
#' @param query_col color in plot
#' @param return_data return data or  plot. defaults to FALSE which means the plot will be printed
#'
#' @return ggplot or data behind plot
#'
#' @export
#'
#' @import Seurat dplyr ggplot2 cowplot stringr
#'
#' @examples
#'
#'

plot_propagation_stats = function(query_seura_object,reference_seurat,label_col_query,label_col,ratio_cut=2,text_size=10,plot_abs = FALSE,n_top = 10,ref_col = "#eda09a",query_col = "#9aa9ed",return_data=FALSE){

  # TODO: add checks ?

  ## create counts per cluster (in ref and query)
  query_count = query_seura_object@meta.data %>% dplyr::group_by(!!rlang::sym(label_col_query)) %>% dplyr::rename(group = !!rlang::sym(label_col_query)) %>% dplyr::count(name="query_count")
  reference_count = reference_seurat@meta.data %>% dplyr::group_by(!!rlang::sym(label_col)) %>% dplyr::rename(group = !!rlang::sym(label_col))  %>% dplyr::count(name="reference_count")
  # combine
  count_both = dplyr::left_join(reference_count,query_count,by = c("group"="group"))
  count_both$query_count[is.na(count_both$query_count)] = 0
  # add some more stats
  count_both = count_both %>% ungroup() %>% dplyr::mutate(reference_pct = round(reference_count / sum(reference_count),5)*100, query_pct = round(query_count / sum(query_count),5)*100) %>%
    dplyr::mutate(relative_log2_ratio = log2(query_pct / reference_pct), absolute_log2_ratio = log2(query_count / reference_count))

  # update dataframe with some general values
  count_both$classification[count_both$relative_log2_ratio >= 0] = "query enriched"
  count_both$classification[count_both$relative_log2_ratio < 0] = "query depleted"
  count_both$relative_log2_ratio[is.infinite(count_both$relative_log2_ratio) & count_both$relative_log2_ratio> 0] = 10
  count_both$relative_log2_ratio[is.infinite(count_both$relative_log2_ratio) & count_both$relative_log2_ratio< 0] = -10

  # plot ratio distribution to show
  maxscales = c(-1*max(abs(count_both$relative_log2_ratio)),max(abs(count_both$relative_log2_ratio)))
  p_ratio = ggplot2::ggplot(count_both,aes(x =relative_log2_ratio,fill=classification)) +
    # geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),fill = ref_col, alpha = alpha_rect) +
    # geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),fill = query_col, alpha = alpha_rect) +
    geom_histogram(bins = 30)+#xlim(maxscales) +
    geom_vline(xintercept = 0) + geom_vline(xintercept = log2(ratio_cut),color="black",linetype="dashed")  +
    geom_vline(xintercept = -1*log2(ratio_cut),color="black",linetype="dashed")+
    xlab("Log2-ratio of cluster pct in reference and query")+ylab("Number of clusters")+
    theme(text = element_text(size=text_size),legend.title = element_blank(),legend.text = element_text(size = text_size/2)) +coord_flip()+scale_fill_manual(values = c("query depleted"=ref_col,"query enriched"=query_col))+
    guides(fill = guide_legend(override.aes = list(size = 0.5)))+
    ggtitle("Ratio of all clusters")
  p_ratio

  ## prepare barplots for clusters
  if(plot_abs){
    use_cols = c("reference_count","query_count")
    ytitle = "Number of cells in cluster"
    col_values = c("reference_count"=ref_col,"query_count"=query_col)
  }else{
    use_cols = c("reference_pct","query_pct")
    ytitle = "Percentage of cells in cluster compared to all cells"
    col_values = c("reference_pct"=ref_col,"query_pct"=query_col)
  }

  ## a barplot with topend different clusters (enriched in query)
  query_enriched_clusters = count_both %>% dplyr::filter(relative_log2_ratio > log2(ratio_cut)) %>% dplyr::slice_max(order_by = relative_log2_ratio,n=n_top) %>%
    tidyr::gather(key = "coltype","value",-group,-classification) #%>5 dplyr::mutate(value = as.numeric(value))
  p1=ggplot(query_enriched_clusters %>% dplyr::filter(coltype %in% use_cols) %>% dplyr::mutate(group = stringr::str_wrap(group,width = 15)),aes(x=group,y=value,group=coltype,fill=coltype))+
    geom_bar(stat="identity", position = "dodge")+coord_flip()+
    ylab(ytitle)+xlab("Clusters")+
    theme(text = element_text(size=text_size))+ guides(fill=guide_legend(title="Dataset"))+
    scale_fill_manual(values = col_values)+
    ggtitle("Top enriched clusters")

  ## a barplot with topend different clusters (enriched in reference)
  reference_enriched_clusters = count_both %>% dplyr::filter(relative_log2_ratio < -1*log2(ratio_cut) & relative_log2_ratio>-10) %>%
    dplyr::slice_min(order_by = relative_log2_ratio,n=n_top)  %>%
    tidyr::gather(key = "coltype","value",-group,-classification)
  p2=ggplot(reference_enriched_clusters %>% dplyr::filter(coltype %in% use_cols) %>% dplyr::mutate(group = stringr::str_wrap(group,width = 15)),aes(x=group,y=value,group=coltype,fill=coltype))+
    geom_bar(stat="identity", position = "dodge")+coord_flip()+
    ylab(ytitle)+xlab("Clusters")+
    theme(text = element_text(size=text_size))+ guides(fill=guide_legend(title="Dataset"))+
    scale_fill_manual(values = col_values)+
    ggtitle("Top depleted clusters")

  if(return_data){
    return(count_both)
  }else{
    #legend_p <- cowplot::get_legend(p_ratio)
    cowplot::plot_grid(p1+ guides(fill=FALSE),p2+ guides(fill=FALSE)+ylab(""),p_ratio,nrow=1)
  }
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
    warning("The networkD3 package must be installed to use this function")
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
