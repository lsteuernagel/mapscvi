
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
### avg_neighbor_distances
##########

#' Estimate quality of mapped data based on distances of between reference cells
#'
#' Calculates average distance between reference neighbors.
#'
#' @param neighbors_object neighbors_object. if NULL tries to find it in query_seurat_object
#' @param reference_map_reduc latent space from reference with same order as idx in `neighbors_object@nn.idx`
#' @param query_seurat_object query seurat. defaults to NULL
#' @param query_nn character vector with name of nn object in query (that contains neighbor indices in reference). only used when neighbors_object is NULL
#' @param distance.metric distance metric, 'cosine' or any available for stats::dist like 'euclidean'
#' @param add_to_seurat add result to query seurat object ? requires query_seurat_object to be  defaults to FALSE
#'
#' @return query_seurat_object with 'avg_neighbor_distance' in meta.data or vector
#'
#' @export
#'
#' @import SeuratObject Seurat dplyr
#'
#' @examples

avg_neighbor_distances = function(neighbors_object=NULL,reference_map_reduc,query_seurat_object=NULL,query_nn = "query_ref_nn",distance.metric="cosine",add_to_seurat=TRUE){

  # get neighbors_object from query seurat if necessary
  if(is.null(neighbors_object)){
  if(query_nn %in%  names(query_seurat_object@neighbors)){
    message("Found ",query_nn)
    neighbors_object = query_seurat_object@neighbors[[query_nn]]
  }else{
    stop("Please provide a valid @neighbors object via the query_nn parameter that labels query neighbors in the reference. Or directly via the neighbors_object parameter")
  }
  }
  # extract nn.idx
  query_to_ref_idx = neighbors_object@nn.idx
  query_to_ref_dists = neighbors_object@nn.dist

  # calculate all pairwise differences of neighbors in reference:
  message("Calculating average distances ...")
  mean_between_distance_reference = apply(query_to_ref_idx,1,function(idx_of_neighbor,latent_space,method =distance.metric){
    if(method == "cosine"){
      #Convert to cosine dissimilarity matrix (distance matrix).
      #https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
      fast_cosine = function(mat){
        sim <- mat / sqrt(rowSums(mat * mat)) # normalize by sum of each row (vector)
        sim <- sim %*% t(sim) #
        dissim <- 1 - sim
        return(dissim)
      }
      distances_of_ref_neighbors = fast_cosine(base::as.matrix(latent_space[idx_of_neighbor,]))
    }else{
      # use R dist for distances: e.g euclidean
      distances_of_ref_neighbors = base::as.matrix(stats::dist( latent_space[idx_of_neighbor,],method = method))
    }
    # get mean and return
    base::mean(distances_of_ref_neighbors[base::upper.tri(distances_of_ref_neighbors)])
  },latent_space=reference_map_reduc@cell.embeddings)
  # return
  if(add_to_seurat & !is.null(query_seurat_object)){
    query_seurat_object@meta.data$avg_neighbor_distance = mean_between_distance_reference
    return(query_seurat_object)
  }else{
    return(mapping_results)

  }
}




#####################
### Older: might remove !?



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

  # get average distance to all neighbors of query cells
  query_dists = apply(query_to_ref_dists[which(query_seurat_object@meta.data$query),2:ncol(query_to_ref_dists)],1,stats::median)
  query_dists = data.frame(Cell_ID = query_seurat_object@meta.data$Cell_ID[query_seurat_object@meta.data$query], median_dist_to_neighbors = query_dists)

  # for each query cell: get average dist of its reference neighbors
  query_neighbor_dists = apply(query_to_ref_idx[which(query_seurat_object@meta.data$query),2:ncol(query_to_ref_idx)],1,function(x,nn_dist){
    stats::median(apply(nn_dist[x,2:ncol(nn_dist)],1,stats::median))
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

