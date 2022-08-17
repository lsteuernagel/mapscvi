# # #####
# # ## load full map
# # #####
#
# library(Seurat)
# library(tidyverse)
#
# #
# hypoMap_v2_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/hypoMap_v2.rds")
#
# # # make new seurat object with reduced size
# reference_hypoMap = CreateSeuratObject(hypoMap_v2_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = hypoMap_v2_seurat@meta.data)
# # add reductions back in
# reference_hypoMap@reductions = hypoMap_v2_seurat@reductions
# # overwrite counts with empty sparse matrix to save space
# empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = reference_hypoMap@assays$RNA@counts@Dim, dimnames = reference_hypoMap@assays$RNA@counts@Dimnames)
# empty_matrix <- as(empty_matrix, "dgCMatrix")
# reference_hypoMap = SetAssayData(object = reference_hypoMap, slot = "counts", new.data = empty_matrix)
# reference_hypoMap = SetAssayData(object = reference_hypoMap, slot = "data", new.data = empty_matrix)
# # overwrite other slots with empty dummy matrix
# dummy=matrix(data = as.numeric())
# reference_hypoMap@assays$RNA@scale.data = dummy
# reference_hypoMap@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
#
# #### downsample
# global_seed = 12345
# cluster_downsample = scUtils::downsample_balanced_iterative(reference_hypoMap@meta.data,target_sample_size = 100000,predictor_var = "C465",stepsize = 5000,global_seed = global_seed)
# table(cluster_downsample$C465)
#
# # subset
# reference_hypoMap_downsample = subset(reference_hypoMap,cells = cluster_downsample$Cell_ID)
#
# # make samller metadata
# reference_hypoMap_downsample_meta = reference_hypoMap_downsample@meta.data %>% dplyr::select(Cell_ID,Dataset,Batch_ID,SRA_ID,nCount_RNA,nFeature_RNA,Technology,Age,Diet,Technology,Strain,
#                                                                              sex = inferred_sex,Region_summarized,Author_Class_Curated, Author_CellType,
#                                                                              C2_named,C7_named,C25_named,C66_named,C185_named,C286_named,C465_named)
# rownames(reference_hypoMap_downsample_meta) = reference_hypoMap_downsample_meta$Cell_ID
# reference_hypoMap_downsample@meta.data = reference_hypoMap_downsample_meta
#
# # Relcalc UMAP:
# reference_hypoMap_downsample = RunUMAP(reference_hypoMap_downsample,
#                                    reduction = "scvi",
#                                    seed.use= 123456,
#                                    dims=1:ncol(reference_hypoMap_downsample@reductions[["scvi"]]@cell.embeddings),
#                                    reduction.name=paste0("umap_","scvi"),
#                                    reduction.key = paste0("umap_","scvi"),
#                                    verbose=TRUE,
#                                    n.neighbors = 25,
#                                    return.model = TRUE)
#
# DimPlot(reference_hypoMap_downsample,group.by = "C465_named")+NoLegend()
#
# # save
# save(reference_hypoMap_downsample,
#      file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_small.RData",
#      compress="xz",
#      compression_level = "9")
#
# # set model
# model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/hypoMap_harmonized_scVI_model/"
# system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))
#
#
#
