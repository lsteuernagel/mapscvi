% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping_functions.R
\name{propagate_labels_prob}
\alias{propagate_labels_prob}
\title{Propagates cell labels base don provided vector using the nearest neighbors in the reference.}
\usage{
propagate_labels_prob(
  neighbors_object = NULL,
  label_vec,
  query_seurat_object = NULL,
  reference_seurat_object = NULL,
  reduction_name_query = "scvi",
  reduction_name_reference = "scvi",
  annoy.metric = "cosine",
  k.param = 30,
  with_euclidean = FALSE,
  add_entropy = FALSE,
  add_to_seurat = FALSE
)
}
\arguments{
\item{neighbors_object}{neighbors-object with reference neighbors of query (see SeuratObject). If NULL will run neighbor detection between query_seurat_object and reference_seurat_object. Requires slot nn.idx and nn.dist.}

\item{label_vec}{vector with labels}

\item{query_seurat_object}{seurat with query data (and a reduction that is projected from reference_seurat_object) to find shared neighbors. NULL by default (not used when neighbors_object is provided)}

\item{reference_seurat_object}{seurat to map onto. NULL by default (not used when neighbors_object is provided)}

\item{reduction_name_query}{name of reduction in query (not used when neighbors_object is provided)}

\item{reduction_name_reference}{name of reduction in reference (not used when neighbors_object is provided)}

\item{annoy.metric}{euclidean or cosine (not used when neighbors_object is provided)}

\item{k.param}{k param for neighbor finding (not used when neighbors_object is provided)}

\item{with_euclidean}{Apply gaussian kernel to smooth distances .only use when distance = euclidean ! (or a neighbors_object based on euclidean distances is provided)}

\item{add_entropy}{re-run to calculate entropy as uncertainty measure}

\item{add_to_seurat}{add to query_seurat_object or return dataframe, requires query_seurat_object to be provided}
}
\value{
seuratobject or dataframe with label propagation results and qc
}
\description{
Cell probabilities similar to scARches algorithm
This function runs per cell
}
