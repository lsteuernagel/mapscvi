
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
#' @param bg_col if non null, the reference will take this color and the query ill be overlayed colored by label_col_query. This will ignore overlay_color and overlay_alpha
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

plot_query_labels = function(query_seura_object,reference_seurat,label_col,label_col_query = "predicted", overlay = FALSE, bg_col = "grey80", overlay_color = "red", overlay_alpha = 0.5,query_pt_size=NULL, query_umap = "umap_scvi",reference_umap="umap_scvi",labelonplot=TRUE,noaxes=TRUE,nolegend=TRUE,...){

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
    # save and remove geom_text layer
    if(labelonplot){
      save_geom_text = p_full$layers[[2]]
      p_full$layers[[2]] =NULL
    }
    # recreate plot if all points are bg col
    if(!is.null(bg_col)){
      reference_seurat$dummy = NA
      p_full=DimPlot(reference_seurat,group.by = "dummy")+scale_color_manual(values = bg_col,na.value=bg_col)+NoLegend()
    }
    if(noaxes){p_full = p_full+Seurat::NoAxes()}
    if(nolegend){ p_full = p_full+Seurat::NoLegend()}
    # adjust alpha
    p_full[[1]]$layers[[1]]$aes_params$alpha = min(overlay_alpha,1)
    # get pt size
    if(is.null(query_pt_size)){
      pt_size = p_full[[1]]$layers[[1]]$aes_params$size
    }else{
      pt_size = query_pt_size
    }
    # plot query points on top
    if(!is.null(bg_col)){
      # if bg color is set we are plotting the color with the query points
      # do the label_col and label_col query overlap ? then use color scale from reference
      if(length(intersect(unique(query_seura_object@meta.data[,label_col_query]),unique(reference_seurat@meta.data[,label_col])))>0){
        test_plot = DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap) # testplot from full data
        testplot_build=ggplot_build(test_plot)$data[1][[1]] # dataframe with colors
        color_mapping_df=as.data.frame(cbind(testplot_build[,"colour"],reference_seurat@meta.data[,c(label_col)])) %>% dplyr::distinct(V1,V2) # make a df with label_col and colours
        color_mapping <- as.character(color_mapping_df$V1) # convert to a named vector for scale_color_manual
        names(color_mapping) <- color_mapping_df$V2
        # add points to plot
        p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2],color=label_col_query),size=pt_size)+
          scale_color_manual(values = color_mapping,na.value= bg_col)
      }else{# if not use default mapping
        p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2],color=label_col_query),size=pt_size)
      }
      p_full = p_full+ggtitle(label_col_query)
    }else{
      # if no bg color use overlay color instead
      p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2]),size=pt_size,color=overlay_color)
    }
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
#' @importFrom tidyr gather
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
    theme(text = element_text(size=text_size),legend.title = element_blank(),legend.text = element_text(size = text_size)) +coord_flip()+scale_fill_manual(values = c("query depleted"=ref_col,"query enriched"=query_col))+
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
  query_enriched_clusters = count_both %>% dplyr::filter(relative_log2_ratio > log2(ratio_cut)) %>% dplyr::arrange(relative_log2_ratio) %>% dplyr::slice_max(order_by = relative_log2_ratio,n=n_top) %>%
    tidyr::gather(key = "coltype","value",-group,-classification) #%>5 dplyr::mutate(value = as.numeric(value))
  p1=ggplot(query_enriched_clusters %>% dplyr::filter(coltype %in% use_cols) %>% dplyr::mutate(group = stringr::str_wrap(group,width = 15)),aes(x=group,y=value,group=coltype,fill=coltype))+
    geom_bar(stat="identity", position = "dodge")+coord_flip()+
    ylab(ytitle)+xlab("Clusters")+
    theme(text = element_text(size=text_size))+ guides(fill=guide_legend(title="Dataset"))+
    scale_fill_manual(values = col_values)+
    ggtitle("Top enriched clusters")

  ## a barplot with topend different clusters (enriched in reference)
  reference_enriched_clusters = count_both %>% dplyr::filter(relative_log2_ratio < -1*log2(ratio_cut) & relative_log2_ratio>-10 )%>% dplyr::arrange(desc(relative_log2_ratio)) %>%
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
### compare_clustering
##########

#' Compare two clustering results
#'
#' For example between propagated and original clusters!
#' TODO: add description
#'
#' @param query_seura_object query seurat object
#' @param clustering_1 column name in query_seura_object metadata with cluster labels
#' @param clustering_2 column name in query_seura_object metadata with cluster labels
#' @param min_cells minimum shared cells to construct a relationship (edge) between two clusters from different clusterings
#' @param min_pct min percentage to construct a relationship (edge) between two clusters from different clusterings
#' @param return_data return data or  plot. defaults to FALSE which means the plot will be printed
#'
#' @return ggplot or data behind plot
#'
#' @export
#'
#' @import Seurat dplyr ggplot2 cowplot stringr
#'
#' @importFrom igraph graph_from_data_frame components degree
#'
#' @examples
#'
#'

# test
# clustering_1 = "predicted_K169_named"
# clustering_2 = "seurat_clusters"

compare_clustering = function(query_seura_object,clustering_1,clustering_2,min_cells = 10,min_pct = 0.1,return_data=FALSE){

  # group by cluster to make summary df
  overview = query_seura_object@meta.data %>% dplyr::select(clustering_1 = !!rlang::sym(clustering_1),clustering_2 = !!rlang::sym(clustering_2)) %>%
    dplyr::group_by(clustering_1) %>% dplyr::add_count(name = "clustering_1_total") %>%
    dplyr::group_by(clustering_2) %>% dplyr::add_count(name = "clustering_2_total") %>%
    dplyr::group_by(clustering_2,clustering_1) %>%
    dplyr::add_count(name="n") %>% dplyr::distinct(clustering_2,clustering_1,.keep_all=TRUE) %>%
    dplyr::mutate(pct_clustering_1 = n / clustering_1_total, pct_clustering_2 = n / clustering_2_total )

  # make a weighted directed graph
  overview_edges = overview %>%  dplyr::select(from = clustering_1, to = clustering_2, weight = pct_clustering_1,n)
  overview_edges = rbind(overview_edges,overview %>% dplyr::select(from = clustering_2, to = clustering_1, weight = pct_clustering_2,n))
  overview_edges =  overview_edges %>% dplyr::filter(n > min_cells, weight > min_pct)
  overview_graph = igraph::graph_from_data_frame(overview_edges,directed = TRUE)
  # find degree
  node_degree = data.frame(node_id = names(igraph::degree(overview_graph, mode = "out")),
                           out_degree = igraph::degree(overview_graph, mode = "out"),
                           in_degree = igraph::degree(overview_graph, mode = "in"),
                           component_id = igraph::components(overview_graph,mode="strong")$membership)

  # add to overview
  overview_relevant = overview %>% dplyr::filter(n > min_cells, pct_clustering_1 > min_pct | pct_clustering_2 > min_pct)
  overview_relevant = dplyr::left_join(overview_relevant,node_degree, by=c("clustering_1"="node_id"))
  overview_relevant = dplyr::left_join(overview_relevant,node_degree[,c("node_id","out_degree","in_degree" )], by=c("clustering_2"="node_id"),suffix = c("_1","_2"))

  # need to find 1:n , m:1 ,m:n
  clustering_1_flagged = overview_relevant %>% ungroup() %>% dplyr::filter(out_degree_1 > 1  & in_degree_1 > 1) %>%
    dplyr::distinct(clustering_1,out_degree_1,in_degree_1) %>% dplyr::arrange(desc(out_degree_1),desc(in_degree_1))

  clustering_2_flagged = overview_relevant %>% ungroup() %>% dplyr::filter(out_degree_2 > 1  & in_degree_2 > 1) %>%
    dplyr::distinct(clustering_2,out_degree_2,in_degree_2) %>% dplyr::arrange(desc(out_degree_2),desc(in_degree_2))

  both_flagged = overview_relevant %>% ungroup() %>% dplyr::filter(out_degree_1 > 1  & in_degree_1 > 1 & out_degree_2 > 1  & in_degree_2 > 1) %>%
    dplyr::distinct(clustering_1,clustering_2,out_degree_1,in_degree_1,out_degree_2,in_degree_2) %>% dplyr::arrange(desc(out_degree_1),desc(in_degree_1),desc(out_degree_2),desc(in_degree_2))

  # return in what form ?
  if(return_data){
    return(overview_relevant)
  }else{
    return(results_flagged = list(clustering_1_flagged = clustering_1_flagged,clustering_2_flagged=clustering_2_flagged,both_flagged=both_flagged))
  }

}

##########
### visNetwork_clustering
##########

#' Compare two clustering results with visNetwork
#'
#' For example between propagated and original clusters!
#' TODO: add description
#'
#' @inheritParams compare_clustering
#' @inheritParams plot_sankey_comparison
#' @param px_height height of interactve visNetwork in pixles
#' @param color_by defaults to NULL. if not NULL a input_clusters column name to color edges by
#' @param palette = Rcolorbrewer palette name
#'
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

visNetwork_clustering = function(query_seura_object,clustering_1,clustering_2,min_cells = 10,min_pct = 0.1,text_size=25,px_height="800px",col1="#cc2118",col2="#302ac9",color_by = NULL,palette = "Set3",return_data=FALSE){

  # optional use of packages:
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    warning("The visNetwork package must be installed to use this function")
    return(NULL)
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    warning("The visNetwork package must be installed to use this function")
    return(NULL)
  }

  # if(length(setdiff(c("clustering_1","clustering_2","pct_clustering_1","pct_clustering_2"),colnames(query_seura_object)))>0){
  #   stop("Cannot find all required column names in query_seura_object: ",paste0(setdiff(c("clustering_1","clustering_2","pct_clustering_1","pct_clustering_2"),colnames(query_seura_object)),collapse=" ; "))
  # }

  # run
  overview = compare_clustering(query_seura_object,clustering_1,clustering_2,min_cells = 10,min_pct = 0.1,return_data=TRUE)
  print(head(overview))

  # make edgelist
  if(is.null(color_by)){color_by=""}
  overview_edges = overview %>%  dplyr::select(from = clustering_1, to = clustering_2, value = pct_clustering_1,n,group = !!rlang::sym(color_by))
  overview_edges = rbind(overview_edges,overview %>% dplyr::select(from = clustering_2, to = clustering_1, value = pct_clustering_2,n,group = !!rlang::sym(color_by)))

  ###plot network
  edges_a = as.data.frame(overview_edges)
  edges_a$arrows = "to"
  if("group" %in% colnames(edges_a)){
    # get palette
    max_col_palette = RColorBrewer::brewer.pal.info[palette,"maxcolors"]
    # creat color mapping df
    col_mapping = data.frame(group_id = unique(edges_a$group))
    col_mapping$color = rep(RColorBrewer::brewer.pal(max_col_palette, "Set3"),ceiling(nrow(col_mapping) / max_col_palette))[1:nrow(col_mapping)] # add idx
    # add to edges
    edges_a = edges_a %>% dplyr::left_join(col_mapping,by=c("group"="group_id"))
  }else{
    edges_a$color = "#878787"
  }

  #nodes
  nodes_a = data.frame(id = unique(c(edges_a$from,edges_a$to)))
  #nodes_a = edges_a %>% dplyr::distinct(from,group) %>% dplyr::select(id = from, group)
  nodes_a$label = nodes_a$id
  nodes_a$font.size = text_size
  nodes_a$color = col1
  nodes_a$color[nodes_a$id %in% overview$clustering_2] = col2

  # render
  visNetwork::visNetwork(nodes_a, edges_a, width = "100%",height = px_height) %>%
    visNetwork::visIgraphLayout(randomSeed = 123,physics = FALSE) %>% visNetwork::visEdges(smooth = TRUE) #%>% visPhysics(stabilization = FALSE)

}


##########
### plot_sankey_comparison
##########

#' Plot comparison between reference and query using sankey plots
#'
#' TODO: add description
#'
#' @param input_clusters edgelist-like dataframe with first two clusters containing ids of clusters
#' @param clustering_1_filter character vector, which clusters to filter on
#' @param clustering_2_filter character vector, which clusters to filter on
#' @param value_col string, column name that specifies values of connection size in sankey
#' @param text_size numeric, font size in sankey
#' @param col1 string, color name for clustering_1
#' @param col2 string, color name for clustering_2
#' @param return_data return list with nodes and edges instead of plotting sankey
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
# input_clusters = overview_relevant
# min_pct = 0.1
# input_clusters = overview %>% dplyr::filter(pct_clustering_1 > min_pct | pct_clustering_2 > min_pct)
# clustering_1_filter = "Slc32a1.Hmx2.Hmx3.Prok2.Th"
# clustering_2_filter = c("57","15","24","29","34")
#
# clustering_1_filter = NULL
# clustering_2_filter = c("40","32","0")

plot_sankey_comparison = function(input_clusters,clustering_1_filter = NULL,clustering_2_filter = NULL,value_col="n",text_size=10,col1="#cc2118",col2="#302ac9",return_data=FALSE){

  # optional use of packages:
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    warning("The networkD3 package must be installed to use this function")
    return(NULL)
  }
  # check that either filter is provided
  if(is.null(clustering_1_filter) & is.null(clustering_2_filter)){
    warning("Please provide a filter for one of the clusterings")
    return(NULL)
  }
  if(colnames(input_clusters)[1] != "clustering_1" | colnames(input_clusters)[2] != "clustering_2"){
    warning("Warning: The first two columns of input_clusters have to be named clustering_1 and clustering_1! Overwriting, but please ensure that they correspond to the cluster ids or names.")
    #return(NULL)
    colnames(input_clusters)[1:2] =c("clustering_1","clustering_2")
  }
  # make edges
  if(!is.null(clustering_1_filter) & !is.null(clustering_2_filter)){
    sankey_edges = input_clusters %>% dplyr::filter(clustering_1 %in% clustering_1_filter | clustering_2 %in% clustering_2_filter)
  }else if(!is.null(clustering_1_filter)){
    sankey_edges = input_clusters %>% dplyr::filter(clustering_1 %in% clustering_1_filter)
  }else{
    sankey_edges = input_clusters %>% dplyr::filter(clustering_2 %in% clustering_2_filter)
  }
  # nodes
  sankey_nodes <- data.frame(name=c(as.character(sankey_edges$clustering_1),as.character(sankey_edges$clustering_2)) %>% unique())
  sankey_nodes$group = "clustering_1"
  sankey_nodes$group[sankey_nodes$name %in% input_clusters$clustering_2] = "clustering_2"
  sankey_edges = as.data.frame(sankey_edges)
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  sankey_edges$IDsource <- match(sankey_edges$clustering_1, sankey_nodes$name)-1
  sankey_edges$IDtarget <- match(sankey_edges$clustering_2, sankey_nodes$name)-1

  # return
  if(return_data){
    return(list(nodes = sankey_nodes, edges = sankey_edges))
  }else{
    #plot
    p <- networkD3::sankeyNetwork(Links = sankey_edges, Nodes = sankey_nodes,
                                  Source = "IDsource", Target = "IDtarget",
                                  Value = value_col, NodeID = "name",NodeGroup = "group",
                                  colourScale = networkD3::JS(paste0('d3.scaleOrdinal() .domain(["clustering_1","clustering_2"]) .range(["',col1,'","',col2,'"]);')),
                                  sinksRight=FALSE,fontSize=text_size)
    p
  }

}
