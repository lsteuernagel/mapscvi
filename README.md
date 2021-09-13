
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mapscvi

<!-- badges: start -->
<!-- badges: end -->

Map single cell expression data in a seurat object into reference scvi
latent space and reference umap based on seurat.

## Installation

Install mapscvi using:

``` r
devtools::install_github("lsteuernagel/mapscvi")
```

In order to use the package python and scvi &gt;= 0.8, as well as R and
Seurat &gt; 4.0.0 are required.

A docker image which comes with a compatible R, Seurat v4, pytorch and
scvi installation can be found here:
<https://hub.docker.com/r/lsteuernagel/r_scvi>

## Example

This package allows to embed new single-cell data stored in Seurat
objects into HypoMap.

The functions used to do this can also be used to embed data into other
models.

An example workflow that emebds the romanov et al. smart-seq dataset
into the HypoMap:

``` r
library(mapscvi)
```

We load the example data which contains a Seurat object.

``` r
load("/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/query_romanov_object.RData")
query_romanov
#> An object of class Seurat 
#> 21143 features across 845 samples within 1 assay 
#> Active assay: RNA (21143 features, 1500 variable features)
```

The test data does not contain any dimensional reductions for
visualization or annotation.

``` r
names(query_romanov@reductions)
#> NULL
```

To map data to the reference HypoMap, we load a reduced Seurat object
with relevant information.

``` r
load("/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/reference_hypoMap.RData")
```

This wrapper function executes all required mapping steps. We provide
the query object as well as a column from HypoMap which we want to
predict.

Take a look at the top clusters that were found in the query:

``` r
head(sort(table(query_romanov@meta.data$predicted),decreasing = TRUE),n = 10)
#> 
#>           Slc17a6.Nrn1.Tbr1.Shox2    Slc17a6.Foxb1.Pitx2.Sepp1.Mobp 
#>                               141                                61 
#>             Slc17a6.Nrn1.Sim1.Trh Slc17a6.Foxb1.Pitx2.Sepp1.Slc7a10 
#>                                53                                49 
#>                         Oxt.Smim3        Slc32a1.Hmx2.Hmx3.Prok2.Th 
#>                                45                                25 
#>        Slc17a6.Nrn1.Sim1.Ebf3.Crh       Slc32a1.Arx.Gad2.Sp9.Fbxw13 
#>                                23                                22 
#>                  Slc17a6.Nrn1.Sst             Slc32a1.Arx.Gad2.Sncg 
#>                                20                                18
```

The package provides plotting functions to visualize the query on the
reference:

Plot them side by side. Here we set labelonplot to False to prevent
cluster labels from being plotted.

``` r
plot_query_labels(query_seura_object=query_romanov,reference_seurat=reference_hypoMap,label_col="K169_named",overlay = FALSE,labelonplot = FALSE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Overlay them query over the reference. The overlay parameters allow to
change the behavior of query points. Here we set labelonplot to True to
include cluster labels on the plot. We can use the Seurat::DimPlot
parameters to further adjust the plots. E.g. by decreasing the size of
the labels.

``` r
plot_query_labels(query_seura_object=query_romanov,reference_seurat=reference_hypoMap,label_col="K169_named",overlay = TRUE,query_pt_size = 0.4,labelonplot = TRUE,label.size=1)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

## Detailed walkthrough

This section contains a more detailed walkthrough of the functions that
are executed when calling ‘map\_new\_seurat\_hypoMap’.

For this we will embed with the mouse LaMannoBrainData stored in a
SingleCellExperiment (available via the scRNAseq package).

``` r
load("/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/sce_lamanno_da.RData")
sce_lamanno_da
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> class: SingleCellExperiment 
#> dim: 18219 243 
#> metadata(0):
#> assays(1): counts
#> rownames(18219): Rp1 Sox17 ... Gm20826_loc2 Erdr1
#> rowData names(0):
#> colnames(243): 1772072122_A04 1772072122_A05 ... 1772099011_H05
#>   1772099012_E04
#> colData names(2): CELL_ID Cell_type
#> reducedDimNames(0):
#> altExpNames(0):
```

The prepare\_query function is able to load Seurat, SingleCellExperiment
or matrix objects for mapping.

Here we are using the prepare\_query\_hypoMap that extends the more
generic prepare\_query function with code to to automatically add
metadata variables that are required for mapping to the HypoMap model.
We also normalize the data, however for the scvi mapping raw counts are
required!

``` r
lamanno_seurat_object = prepare_query_hypoMap(sce_lamanno_da,suffix="lamanno_da",normalize=TRUE)
lamanno_seurat_object
#> An object of class Seurat 
#> 18219 features across 243 samples within 1 assay 
#> Active assay: RNA (18219 features, 0 variable features)
```

This new seurat object is compatible with the downstream functions for
mapping the data.

Next predict\_query can be used to embed the query data into the latent
space of scvi. We have to specify a query path, which defaults to the
HypoMap reference that is part of the package and the number of epochs
for training during mapping (10-20 should be sufficient. TODO: test
this!).

``` r
names(lamanno_seurat_object@reductions)
#> [1] "scvi"
```

The scvi reduction is a pca-like low dimensional space that can be used
to embed the data into the same UMAP as the reference object. This
requires an existsing UMAP model in the reference Seurat object that was
calculated based on teh same scvi latent space.

We load the reference seurat object.

``` r
load("/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/reference_hypoMap.RData")
```

Then we can calculate nearest neighbors and UMAP based on the reference.
Additionally we can project labels (any categorical metadata column from
the reference) using the same nearest neighbor graph. This helps with
consistency between label propagation and UMAP. However there might be
more accurate ways to propagate labels using other classifiers such as a
random forest or scANVI.

To propagate labels with the project\_query function we can provide a
vector of the same lengths as the reference cells (and same order!).
preferabbly a column from the metadata of the reference seurat object.

This can then be used to plot the results:

``` r
head(sort(table(lamanno_seurat_object@meta.data$predicted),decreasing = TRUE),n = 10)
#> 
#>         Slc17a6.Nrn1.Tbr1.Ebf3 Slc17a6.Foxb1.Pitx2.Sepp1.Mobp 
#>                            200                             25 
#>   Slc32a1.Satb2.Fam159b.Slc6a3   Slc17a6.Nrn1.Tbr1.Pitx2.Irx3 
#>                             10                              6 
#>          Slc32a1.Arx.Gad2.Sncg         Slc32a1.Hmx2.Lef1.Prph 
#>                              1                              1
```

``` r
plot_query_labels(query_seura_object=lamanno_seurat_object,reference_seurat=reference_hypoMap,label_col="K169_named",overlay = FALSE,labelonplot = FALSE)
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

``` r
plot_query_labels(query_seura_object=lamanno_seurat_object,reference_seurat=reference_hypoMap,label_col="K169_named",overlay = TRUE,query_pt_size = 0.4,labelonplot = TRUE,label.size=1)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

This obviously didn’t work out very well.

TODO: Add some indicators for quality of mapping.
