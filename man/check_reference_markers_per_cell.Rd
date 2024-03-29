% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_functions.R
\name{check_reference_markers_per_cell}
\alias{check_reference_markers_per_cell}
\title{Estimate quality of mapped data based on marker genes from reference}
\usage{
check_reference_markers_per_cell(
  query_seurat_object,
  marker_genes,
  assay = "RNA",
  query_label_col = "predicted",
  min_expr = 0,
  global_seed = 12345
)
}
\arguments{
\item{query_seurat_object}{Seurat object}

\item{marker_genes}{a named list of character vectors with (marker) gene names and names corresponding to values in query_label_col}

\item{assay}{assay name. defaults to RNA}

\item{query_label_col}{string, column-name in metadata of query_seurat_object. values of have to correspond to names of marker_genes}

\item{min_expr}{min expression to be considered expressed (defaults to 0)}

\item{global_seed}{seed}
}
\value{
query_seurat_object with quality results in metadata
marker_pct: how many markers of projected reference cluster are expressed in each cell
}
\description{
Check marker gene expression in query cells based on predicted clusters and markers from reference.
}
