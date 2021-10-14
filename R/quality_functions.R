##########
### check_distance_neighbors
##########

#' Estimate quality of mapped data based on distances of query cells from reference cells
#'
#' Calculates average distance to reference neighbors and normalized average distances based on average distances of these reference neighbors within the reference dataset.
#'
#' @inheritParams project_query
#' @param reference_seurat reference seurat
#' @param reduction_name name of scvi reduction in reference AND query. Also expect a umap named as 'umap_'reference_reduction
#' @param query_nn character vector with name
#' @param reference_nn character vector with name
#' @param k_param k neighbors
#' @param annoy.metric string passed to Seurat's FindNeighbors
#'
#' @return query_seurat_object with quality results in metadata
#' nn_distance_reference: what is the median distance of a query cell to all neighbor cells in reference
#' normalized_nn_distance: nn_distance_reference / median_dist_refNeighbors (how comparable are the query cells distances to their neighboring cell distances)
#'
#' @export
#'
#' @import SeuratObject Seurat dplyr
#'
#' @examples

check_distance_neighbors = function(query_seurat_object,reference_seurat,reduction_name="scvi",query_nn = "query_ref_nn",reference_nn = "reference_nn",k_param=30,annoy.metric="cosine",global_seed=12345){

  # check that reduction_name exists
  #TODO

  # add query boolean:
  #query_seurat_object@meta.data$query =TRUE

  # query
  if(query_nn %in%  names(query_seurat_object@neighbors)){
    message("Found ",query_nn)
  }else{
    stop("Please provide a valid @neighbors object via the query_nn parameter that labels query neighbors in the reference.")
    #query_seurat_object = Seurat::FindNeighbors(query_seurat_object,reduction =  reduction_name,k.param = k_param,graph.name = query_nn,annoy.metric = annoy.metric,
    #                                    return.neighbor = TRUE,dims = 1:ncol(query_seurat_object@reductions[[reduction_name]]@cell.embeddings))
  }
  query_to_ref_idx = query_seurat_object@neighbors[[query_nn]]@nn.idx
  query_to_ref_dists = query_seurat_object@neighbors[[query_nn]]@nn.dist

  # ref
  if(reference_nn %in%  names(reference_seurat@neighbors)){
    message("Found ",reference_nn)
  }else{
    reference_seurat = Seurat::FindNeighbors(reference_seurat,reduction =  reduction_name,k.param = k_param,graph.name = reference_nn,annoy.metric = annoy.metric,
                                             return.neighbor = TRUE,dims = 1:ncol(reference_seurat@reductions[[reduction_name]]@cell.embeddings))
  }
  ref_idx = reference_seurat@neighbors[[reference_nn]]@nn.idx
  ref_dists = reference_seurat@neighbors[[reference_nn]]@nn.dist
#
#   message("Calculating average distances ...")
#   mean_between_distance_reference = apply(neighbors_reference_object@nn.idx,1,function(idx_of_neighbor,latent_space){
#     distances_of_ref_neighbors = as.matrix(dist( latent_space[idx_of_neighbor,],method = "euclidean"))
#     mean(distances_of_ref_neighbors[upper.tri(distances_of_ref_neighbors)])
#   },latent_space=latent_space)
#
#
#
  # get average distance to all neighbors of query cells
  query_dists = apply(query_to_ref_dists[which(query_seurat_object@meta.data$query),2:ncol(query_to_ref_dists)],1,stats::median)
  query_dists = data.frame(Cell_ID = query_seurat_object@meta.data$Cell_ID[query_seurat_object@meta.data$query], median_dist_to_neighbors = query_dists)

  # for each query cell: get average dist of its reference neighbors
  query_neighbor_dists = apply(query_to_ref_idx[which(query_seurat_object@meta.data$query),2:ncol(query_to_ref_idx)],1,function(x,nn_dist){
    median(apply(nn_dist[x,2:ncol(nn_dist)],1,stats::median))
  },nn_dist=ref_dists)
  # summarise
  query_neighbor_dists = data.frame(Cell_ID = query_seurat_object@meta.data$Cell_ID[query_seurat_object@meta.data$query], median_dist_refNeighbors = query_neighbor_dists)

  # compute normalized score
  query_seurat_object@meta.data$nn_distance_reference = query_dists$median_dist_to_neighbors
  query_seurat_object@meta.data$normalized_nn_distance = query_dists$median_dist_to_neighbors / query_neighbor_dists$median_dist_refNeighbors
  message("Adding result to metadata as: ","nn_distance_reference, ","normalized_nn_distance")

  # return
  return(query_seurat_object)

}

##########
### check_freq_neighbors
##########

#' Estimate quality of mapped data based on neighborhood of query cells
#'
#' Calculates a shared KNN to see how well query cells mix with reference cells. Depends on the number of query cells (many query cells will likely also lead to more query dataset neighbors)
#'
#' @inheritParams project_query
#' @inheritParams check_distance_neighbors
#' @param reference_seurat reference seurat
#'
#' @return query_seurat_object with quality results in metadata
#' query_neighbor_pct: how many neighbors of a query cells come from query dataset. Can either indicate that query contains many more cells than reference or that celltype is absent in reference.
#'
#' @export
#'
#' @import SeuratObject Seurat dplyr
#'
#' @examples

check_freq_neighbors = function(query_seurat_object,reference_seurat,reduction_name="scvi",k_param=30,annoy.metric="cosine",global_seed=12345){

  # check that reduction_name exists
  #TODO

  # add query boolean:
  query_seurat_object@meta.data$query =TRUE
  # merge objects
  merged_object = merge(reference_seurat,query_seurat_object,merge.dr = reduction_name)
  # set query variable to false
  merged_object@meta.data$query[is.na(merged_object@meta.data$query)] =FALSE
  # need to run FindNeighbors
  message("Calculating shared neighbors ...")
  merged_object = Seurat::FindNeighbors(merged_object,reduction =  reduction_name,k.param = k_param,graph.name = "new_nn",annoy.metric = annoy.metric,
                                        return.neighbor = TRUE,dims = 1:ncol(merged_object@reductions[[reduction_name]]@cell.embeddings))
  nn_idx = merged_object@neighbors$new_nn@nn.idx
  nn_dist = merged_object@neighbors$new_nn@nn.dist

  # use nn_idx to get number of query neighbors for each query cells
  message("Calculating neighbor frequencies ...")
  query_freq = apply(nn_idx[which(merged_object@meta.data$query),2:ncol(nn_idx)],1,function(row,query){
    freq_batch = sum(query[row]) / length(row)
  }, query = merged_object@meta.data$query)
  query_freq = data.frame(Cell_ID = merged_object@meta.data$Cell_ID[merged_object@meta.data$query], query_neighbor_pct = query_freq)

  ## add result
  query_seurat_object@meta.data$query_neighbor_freq = query_freq$query_neighbor_pct
  message("Adding result to metadata as: ","query_neighbor_freq")
  # return
  return(query_seurat_object)

}

##########
### check_reference_markers_per_cell
##########

#' Estimate quality of mapped data based on marker genes from reference
#'
#' Check marker gene expression in query cells based on predicted clusters and markers from reference.
#'
#' @inheritParams project_query
#' @inheritParams check_distance_neighbors
#' @param marker_genes a named list of character vectors with (marker) gene names and names corresponding to values in query_label_col
#' @param query_label_col string, column-name in metadata of query_seurat_object. values of have to correspond to names of marker_genes
#' @param min_expr min expression to be considered expressed (defaults to 0)
#'
#' @return query_seurat_object with quality results in metadata
#' marker_pct: how many markers of projected reference cluster are expressed in each cell
#'
#' @export
#'
#' @import SeuratObject Seurat dplyr Matrix
#'
#' @examples

check_reference_markers_per_cell = function(query_seurat_object,marker_genes,assay="RNA",query_label_col="predicted",min_expr=0,global_seed=12345){

  # check input arguments
  if(!query_label_col %in% colnames(query_seurat_object@meta.data)){
    stop("Error: Cannot find '",query_label_col,"' in query_seurat_object@meta.data. Stopping.")
  }else{
    all_labels = unique(query_seurat_object@meta.data[,query_label_col])
  }

  if(length(setdiff(all_labels,names(marker_genes)))>0){
    stop("Warning: Some labels are not part of marker_genes list (cannot find in list names!)")
  }
  # run pct of expressed markers per cell
  res_marker_freq_list = sapply(all_labels,function(label,all_labels,count_matrix,marker_genes){
    label_idx = which(all_labels == label)
    label_genes = marker_genes[[label]]
    if(length(label_genes)>1){
      subset_matrix = count_matrix[rownames(count_matrix) %in% label_genes,label_idx]
      subset_matrix[subset_matrix>min_expr] = 1
      if(length(label_idx)>1){
        freq_df = data.frame(Cell_ID = colnames(subset_matrix), marker_pct = Matrix::colSums(subset_matrix) / length(label_genes))
      }else if(length(label_idx)==1){
        freq_df = data.frame(Cell_ID = colnames(count_matrix)[label_idx], marker_pct = sum(subset_matrix) / length(label_genes))
      }else{
        freq_df =NULL
      }
    }else{
      freq_df =NULL
    }
    freq_df
  },all_labels=query_seurat_object@meta.data[,query_label_col],count_matrix=query_seurat_object@assays[[assay]]@data,marker_genes=marker_genes,simplify = FALSE)
  res_marker_freq = do.call(rbind,res_marker_freq_list)

  ## add to metadata
  if("marker_pct" %in% colnames(query_seurat_object@meta.data)){
    query_seurat_object@meta.data = query_seurat_object@meta.data %>% dplyr::select(-marker_pct)
    rownames(query_seurat_object@meta.data) = query_seurat_object@meta.data$Cell_ID
  }
  query_seurat_object@meta.data = dplyr::left_join(query_seurat_object@meta.data,res_marker_freq,by="Cell_ID")
  rownames(query_seurat_object@meta.data) = query_seurat_object@meta.data$Cell_ID
  message("Adding result to metadata as: ","marker_pct")

  return(query_seurat_object)

}


##########
### adjusted_cell_probabilities
##########

#' Estimate quality of mapped data based on marker genes from reference
#'
#' Cell probabilities similar to scARches algorithm
#' This function runs per cell
#'
#' @param dist_Nc vector of length k with euclidean distances to neigbors
#' @param label_of_neighbor vector of length k with labels of neighbors
#' @param result_type return type: "all","label","entropy"
#'
#' @return label probability depending on result_type
#'
#' @export

adjusted_cell_probabilities = function(dist_Nc,labels_of_neighbor,result_type="all"){

  # step 1: Distances and input
  # ... input of function
  local_labels = unique(labels_of_neighbor)
  k = length(labels_of_neighbor)
  # step 2: compute the standard deviation of the nearest distances
  sd_nc = sqrt(sum(dist_Nc^2) / k )
  # step 3: apply Gaussian kernel to distances
  d_app = exp(-1*(dist_Nc/(2/sd_nc)^2) )
  # step 4: we computed the probability of assigning each label y to the query cell c by normalizing across all adjusted distances
  label_probabilities = tapply(d_app,INDEX = labels_of_neighbor,FUN = sum) / sum(d_app)
  # return result
  if(result_type == "all"){
    return(label_probabilities)
    # or summarise further before returning:
  }else{
    if(result_type == "label"){
      # step 5: uncertainty score: 1 - prob_label --> I rather return directly the highest label with its name
      return(label_probabilities[label_probabilities == max(label_probabilities)])
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

#' Estimate quality of mapped data based on marker genes from reference
#'
#' Cell probabilities similar to scARches algorithm
#' This function runs per cell
#'
#' @param query_seurat_object
#' @param reference_seurat_object v
#' @param label_col
#' @param reduction_name_query
#' @param reduction_name_reference
#' @param k.param
#' @param add_entropy
#' @param add_to_seurat
#'
#' @return object
#'
#' @export

propagate_labels_prob = function(query_seurat_object,reference_seurat_object,label_col,reduction_name_query="scvi",reduction_name_reference="scvi",k.param=30, add_entropy =FALSE, add_to_seurat =TRUE){

  # need euclidean distances neighbors
  neighbors_object = Seurat::FindNeighbors(reference_seurat_object@reductions[[reduction_name_reference]]@cell.embeddings,
                                           query = query_seurat_object@reductions[[reduction_name_query]]@cell.embeddings,
                                           k.param = k.param, return.neighbor =TRUE,
                                           annoy.metric="euclidean")

  # define label vector
  all_labels = reference_seurat_object@meta.data[,label_col]

  # apply max prob per cell function
  message("Estimate probabilities")
  n=nrow(neighbors_object@nn.dist)
  max_probabilities = sapply(1:n,function(x,distances,neighbor_idxs,labels){
    dist_Nc = distances[x,]
    label_of_neighbor = labels[neighbor_idxs[x,]]
    prob = adjusted_cell_probabilities(dist_Nc = dist_Nc,labels_of_neighbor = label_of_neighbor,result_type = "label")
    prob
  },distances = neighbors_object@nn.dist,neighbor_idxs = neighbors_object@nn.idx,labels = all_labels)

  # mapping results
  mapping_results = data.frame(Cell_ID = neighbors_object@cell.names[1:n], predicted = names(max_probabilities), prediction_probability = as.numeric(max_probabilities))

  # not super efficient because I apply the function twice ....
  if(add_entropy){
    message("Estimate entropy")
    # apply max prob per cell function
    mapping_results$prediction_entropy = sapply(1:n,function(x,distances,neighbor_idxs,labels){
      dist_Nc = distances[x,]
      label_of_neighbor = labels[neighbor_idxs[x,]]
      prob = adjusted_cell_probabilities(dist_Nc = dist_Nc,labels_of_neighbor = label_of_neighbor,result_type = "entropy")
      prob
    },distances = neighbors_object@nn.dist,neighbor_idxs = neighbors_object@nn.idx,labels = all_labels)
  }
  # return
  if(add_to_seurat){
    # remove if existing for clean join
    keep_names = colnames(query_seurat_object@meta.data)[!colnames(query_seurat_object@meta.data) %in% c("predicted","prediction_probability","prediction_entropy")]
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


# query_snseq_neurons = propagate_labels_prob(query_seurat_object = query_snseq_neurons,
#                                             reference_seurat_object = neuron_map_seurat,
#                                             label_col = "K329_pruned",k.param = 30,add_to_seurat = TRUE,add_entropy = TRUE)
# query_snseq_neurons$prediction_entropy = 1 - query_snseq_neurons$prediction_entropy
# DimPlot(query_snseq_neurons,group.by = "predicted",reduction = "umap_scvi")+NoLegend()+NoAxes()
# FeaturePlot(query_snseq_neurons,features = "prediction_probability")
# FeaturePlot(query_snseq_neurons,features = "prediction_entropy")

