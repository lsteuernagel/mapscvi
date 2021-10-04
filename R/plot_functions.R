
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
#' @param edgelist TODO
#' @param heatmap_matrix a matrix with rownames
#' @param leaf_level which level to use as leaves ?
#' @param anno_df provide a dataframe which contains mappings of all cluster ids and names in the cluster tree. Defaults to NULL which means automatic construction.
#' @param metadata metadata with cluster ids and names, if anno_df is NULL this has to be provided to construct anno_df
#' @param level_pattern regex for cluster level, if anno_df is NULL this has to be provided to construct anno_df. defaults to 'K[0-9]+'
#' @param cluster_id_pattern regex for cluster id column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_pruned'
#' @param cluster_name_pattern regex for cluster name column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_named'
#' @param label_size label size on tree. defaults to 2
#'
#' @return ggtree object or plot
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#'

plot_cluster_tree = function(edgelist,heatmap_matrix=NULL,leaf_level=NULL,anno_df=NULL,metadata=NULL,level_pattern = "K[0-9]+",cluster_id_pattern = "_pruned",
                             cluster_name_pattern = "_named",label_size = 2){

  # load for tests:
  # edgelist =  neuron_map_seurat@misc$pruned_edgelist
  # leaf_level = 6
  # metadata = neuron_map_seurat@meta.data
  # level_pattern = "K[0-9]+"
  # cluster_id_pattern = "_pruned"
  # cluster_name_pattern = "_named"
  # anno_df=NULL
  heatmap_data = metadata %>% dplyr::select(Cell_ID,K169_named) %>% dplyr::group_by(K169_named) %>%  #dplyr::filter(predicted_Campbell!="NA")
    dplyr::add_count(name = "presence") %>% dplyr::distinct(K169_named,.keep_all=TRUE) %>%dplyr::ungroup() %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup()# %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  # tidyr::spread(key = 1,value=presence)
  heatmap_data = heatmap_data %>% dplyr::full_join(anno_df[,c("cluster_id","cluster_name")],by=c("K169_named"="cluster_name"))  %>% dplyr::filter(grepl("K169",cluster_id))
  heatmap_matrix = as.matrix(heatmap_data[,"presence"])
  rownames(heatmap_matrix) = heatmap_data$cluster_id
  heatmap_matrix[is.na(heatmap_matrix)] = 0

  # optional use of packages: ggtree imports also tidytree and treeio!
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    warning("The ggtree package must be installed to use this function")
    return(NULL)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("The igraph package must be installed to use this function")
    return(NULL)
  }
  # check that req columns exist
  if(length(setdiff(c("from","to","level"),colnames(edgelist)))>0){
    warning("Error: Wrong edgelist format. Requires columns: from, to, level")
    return(NULL)
  }
  # check that leaf_level is there
  if(!leaf_level %in% edgelist$level){
    warning("Error: leaf_level '",leaf_level,"' cannot be found in level column of edgelist")
    return(NULL)
  }

  # construct a dataframe with the required annotations
  if(is.null(anno_df)){
    if(is.null(metadata)){
      warning("Error: Please provide metadata with corresponding cluster ids and names that match provided edgelist")
      return(NULL)
    }
    if(!any(grepl(cluster_id_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_id_pattern '",cluster_id_pattern,"' in provided metadata")
      return(NULL)
    }
    if(!any(grepl(cluster_name_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_name_pattern '",cluster_name_pattern,"' in provided metadata")
      return(NULL)
    }
    # anno_df: cluster id, clean_names, clusterlevel, ncells, first_cluster_name
    pruned_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_id_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = str_extract(colname,level_pattern)) %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    named_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_name_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = str_extract(colname,level_pattern))  %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    both_map = dplyr::left_join(pruned_ids,named_ids,by=c("Cell_ID"="Cell_ID","clusterlevel"="clusterlevel")) %>% dplyr::select(cluster_id = id.x,cluster_name = id.y)

    anno_df = neuron_map_seurat@misc$pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
    anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][length(strsplit(x,"\\.")[[1]])]})
  }else{
    # check that provided anno_df is valid:
    if(length(setdiff(c("cluster_id","clusterlevel","cluster_name","first_cluster_name"),colnames(anno_df)))>0){
      stop("Wrong anno_df format. Required columns: cluster_id, clusterlevel, cluster_name, first_cluster_name")
    }
  }

  # if a heatmap matrix is provided, this function tries to infer the leaflevel based on the matrix
  if(!is.null(heatmap_matrix) & is.null(leaf_level)){
   # TODO
  }

  # reduce edgelist to certain level and from and to cols
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level<=as.numeric(leaf_level),1:2]

  ## convert to treedata
  # only take
  tree_data_igraph = igraph::graph_from_edgelist(as.matrix(edgelist))
  tree_data_phylo = treeio::as.phylo(tree_data_igraph)
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)

  # add labels from annotation_df
  tree_data_tibble = dplyr::left_join(tree_data_tibble,anno_df,by=c("label"="cluster_id"))

  # update additional columns
  tree_data_tibble$first_cluster_name[is.na(tree_data_tibble$first_cluster_name)]=""
  tree_data_tibble$nodesize = 1 # default node size
  tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edgelist) # count children number
  tree_data_tibble$n_siblings = sapply(tree_data_tibble$label,function(x,el){ # count siblings
    parent = el$from[el$to==x]
    return(length(el$to[el$from==parent])-1)
  },el=edgelist)
  tree_data_tibble$tip.label =NA
  tree_data_tibble$tip.label[tree_data_tibble$n_children==0] = tree_data_tibble$node[tree_data_tibble$n_children==0] # add tip labels if leaf
  tree_data_tibble$first_cluster_name[ tree_data_tibble$n_children<2 & is.na(tree_data_tibble$tip.label)] = "" # if only one child node: changeto ""

  # convert back to treedata
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))

  ## heatmap
  # heatmap_data = metadata %>% dplyr::select(Cell_ID,K169_named) %>% dplyr::group_by(K169_named) %>%  #dplyr::filter(predicted_Campbell!="NA")
  #   dplyr::add_count(name = "presence") %>% dplyr::distinct(K169_named,.keep_all=TRUE) %>%dplyr::ungroup() %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup()# %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  # # tidyr::spread(key = 1,value=presence)
  # heatmap_data = heatmap_data %>% dplyr::full_join(anno_df[,c("cluster_id","cluster_name")],by=c("predicted_K169_named"="cluster_name"))  %>% dplyr::filter(grepl("K169",cluster_id))
  # heatmap_matrix = as.matrix(heatmap_data[,"presence"])
  # rownames(heatmap_matrix) = heatmap_data$cluster_id
  # heatmap_matrix[is.na(heatmap_matrix)] = 0

  #plot circular tree
  circular_tree =ggtree::ggtree(tree_data, layout = 'circular', branch.length='none')+
    # geom_text(aes(label=first_cluster_name))
    ggtree::geom_nodelab(aes(x=branch, label=first_cluster_name), size=label_size,vjust=-.5, color="darkred")+
    ggtree::geom_tiplab(aes(x=branch, label=first_cluster_name), size=label_size,vjust=-.5,color="darkred")+
    ggtree::geom_nodepoint(aes(subset = n_children > 1))#+geom_tippoint()

  if(!is.null(heatmap_matrix)){
    heatmap_matrix[is.na(heatmap_matrix)] = 0
    #circular_tree
    circular_tree_heat = ggtree::gheatmap(circular_tree,data=heatmap_matrix,offset = 0.2, width = 0.2,colnames_angle=90)+
      ggplot2::scale_fill_gradient(low = "white",high="darkred",limits=c(0,5),oob=squish) +
      ggplot2::guides(fill=guide_legend(title="Legend"))
    #show:
    circular_tree_heat
  }else{
    circular_tree
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
