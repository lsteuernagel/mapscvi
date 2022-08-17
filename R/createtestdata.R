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
# # remove some metadata to make the object smaller
# #temp_meta = reference_hypoMap@meta.data
# # temp_meta = temp_meta %>% dplyr::select(Cell_ID,Dataset,Batch_ID,Technology,inferred_sex,nFeature_RNA,nCount_RNA,
# #                                         Age,Diet,Author_Class_Curated,Region_summarized,Author_CellType,
# #                                         C2_named,C7_named,C25_named,C66_named,C185_named,C286_named,C465_named)
#
# reference_information = reference_hypoMap@meta.data %>% dplyr::distinct(C465,C465_named,C2_named,C7_named,C25_named,C66_named,C185_named,C286_named,Region_summarized)
# temp_meta = reference_hypoMap@meta.data %>% dplyr::select(Cell_ID,Batch_ID,C465)
# rownames(temp_meta) = temp_meta$Cell_ID
# #reference_hypoMap@meta.data = temp_meta
#
# reference_list_hypomap = list(
#   reference_reduction = reference_hypoMap@reductions$scvi,
#   reference_umap = reference_hypoMap@reductions$umap_scvi,
#   reference_labels = temp_meta,
#   reference_information
# )
#
# # ref info
# save(reference_list_hypomap,
#      file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_list_hypomap.RData",
#      compress="xz",
#      compression_level = "9")
#
# # save
# save(reference_hypoMap,
#      file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap.RData",
#      compress="xz",
#      compression_level = "9")
#
# # set model
# model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap"
# system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))
#
# load("data/reference_hypoMap.RData")

# copie to /beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_v2_objects

# # #####
# # ## load neuron map
# # #####
# map_name = "hypothalamus_neurons_reference"
# map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/"
# map_seurat_path = paste0(map_path,map_name,".rds")
# neuron_map_seurat = readRDS(map_seurat_path)
#
# # # make new seurat object with reduced size
# reference_hypoMap = CreateSeuratObject(neuron_map_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = neuron_map_seurat@meta.data)
# # add reductions back in
# reference_hypoMap@reductions = neuron_map_seurat@reductions
# # overwrite counts with empty sparse matrix to save space
# empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = reference_hypoMap@assays$RNA@counts@Dim, dimnames = reference_hypoMap@assays$RNA@counts@Dimnames)
# empty_matrix <- as(empty_matrix, "dgCMatrix")
# reference_hypoMap = SetAssayData(object = reference_hypoMap, slot = "counts", new.data = empty_matrix)
# reference_hypoMap = SetAssayData(object = reference_hypoMap, slot = "data", new.data = empty_matrix)
# # overwrite other slots with empty dummy matrix
# dummy=matrix(data = as.numeric())
# reference_hypoMap@assays$RNA@scale.data=dummy
# reference_hypoMap@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
# # remove some metadata to make the object smaller
# temp_meta = reference_hypoMap@meta.data
# temp_meta = temp_meta %>% dplyr::select(Cell_ID,Dataset,Batch_ID,Technology,inferred_sex.x,nFeature_RNA,nCount_RNA,rpl_signature_expr_median,
#                                         Age,Diet,Subregion,suggested_region_curated,Author_CellType,K4_named,K14_named,K31_named,K98_named,K169_named,K240_named,
#                                         K329_named,K4_pruned,K14_pruned,K31_pruned,K98_pruned,K169_pruned,K240_pruned,K329_pruned)
# rownames(temp_meta) = temp_meta$Cell_ID
# reference_hypoMap@meta.data = temp_meta
#
# # add markers as a list for K169_pruned
# markers_K169_pruned_list = neuron_map_seurat@misc$markers_comparisons_all %>%
#   dplyr::filter(grepl("K169",cluster_1) & specificity > 0.5 & p_val_adj < 0.001) %>% dplyr::arrange(desc(specificity))
# markers_K169_pruned_list = base::split(markers_K169_pruned_list$gene,markers_K169_pruned_list$cluster_1)
# reference_hypoMap@misc$markers_K169_pruned = markers_K169_pruned_list
#
# # print size and change name
# print(object.size(reference_hypoMap) / 1000000)
# reference_hypoMap_neurons = reference_hypoMap
#
# save(reference_hypoMap_neurons,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_neurons.RData",compress="xz",compression_level = "9")

# # set model
# system(paste0("mkdir -p /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))
# model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_models/hypothalamus_neurons_reference_model/"
# system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))
#
# #####
# ## load full map
# #####
#
# project_name = "hypothalamus_full_map"
# project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/"
# seurat_file_name = paste0(project_path,project_name,".h5Seurat")
# full_map_seurat = SeuratDisk::LoadH5Seurat(seurat_file_name)
#
# # # make new seurat object with reduced size
# reference_hypoMap_full = CreateSeuratObject(full_map_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = full_map_seurat@meta.data)
# # add reductions back in
# reference_hypoMap_full@reductions = full_map_seurat@reductions
# # overwrite counts with empty sparse matrix to save space
# empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = reference_hypoMap_full@assays$RNA@counts@Dim, dimnames = reference_hypoMap_full@assays$RNA@counts@Dimnames)
# empty_matrix <- as(empty_matrix, "dgCMatrix")
# reference_hypoMap_full = SetAssayData(object = reference_hypoMap_full, slot = "counts", new.data = empty_matrix)
# reference_hypoMap_full = SetAssayData(object = reference_hypoMap_full, slot = "data", new.data = empty_matrix)
# # overwrite other slots with empty dummy matrix
# dummy=matrix(data = as.numeric())
# reference_hypoMap_full@assays$RNA@scale.data = dummy
# reference_hypoMap_full@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
#
# # remove some metadata to make the object smaller
# temp_meta = reference_hypoMap_full@meta.data
# temp_meta = temp_meta %>% dplyr::select(Cell_ID,Dataset,Batch_ID,Technology,inferred_sex,nFeature_RNA,nCount_RNA,rpl_signature_expr_median,
#                                         Age,Diet,Subregion,Curated_Class,Curated_CellType,Author_CellType)
# rownames(temp_meta) = temp_meta$Cell_ID
# reference_hypoMap_full@meta.data = temp_meta
#
# # save
# save(reference_hypoMap_full,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_full.RData",compress="xz",compression_level = "9")
#
# # set model
# model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_models/scVI_hypothalamus_full_map_model/"
# system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))
#
#
# #####
# ## load example data
# #####
#
# ## load romanov query
# suffix ="mapped_data_yeo_romanov" # a name
# query_romanov_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/mapped_data_Romanov_neurons/mapped_data_Romanov_neurons.h5Seurat" # seurat object to load
# query_romanov = SeuratDisk::LoadH5Seurat(query_romanov_path)
# query_romanov@reductions = list()
# query_romanov@meta.data = query_romanov@meta.data[,1:20]
# # overwrite other slots with empty dummy matrix
# empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = query_romanov@assays$RNA@counts@Dim, dimnames = query_romanov@assays$RNA@counts@Dimnames)
# empty_matrix <- as(empty_matrix, "dgCMatrix")
# query_romanov = SetAssayData(object = query_romanov, slot = "data", new.data = empty_matrix)
# dummy=matrix(data = as.numeric())
# query_romanov[["RNA"]]@scale.data = dummy
# print(object.size(query_romanov) / 1000000)
#
# # remove some metadata to make the object smaller
# temp_meta = query_romanov@meta.data
# temp_meta = temp_meta %>% dplyr::select(Cell_ID,Sample_ID,nFeature_RNA,nCount_RNA,percent.mt,Sex,Author_CellType,Subregion,predicted_Curated_Class)
# rownames(temp_meta) = temp_meta$Cell_ID
# query_romanov@meta.data = temp_meta
#
# #
# require(scRNAseq)
# sce_lamanno_da <- LaMannoBrainData(which = "mouse-adult",ensembl=FALSE)
# object.size(sce_lamanno_da) / 1000000
# lamanno_seurat = prepare_query(sce_lamanno_da,suffix="lamanno_da",normalize=TRUE)
# lamanno_seurat
# save(lamanno_seurat,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/seurat_lamanno_da.RData")
#
# #save testdata
# save(query_romanov,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/query_romanov.RData",compress="xz",compression_level = "9")
# save(sce_lamanno_da,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/sce_lamanno_da.RData",compress="xz",compression_level = "9")
#
# #####
# ## function for downsampling
# #####
#
# downsample_balanced_iterative <- function(metadata,target_sample_size,predictor_var,stepsize =100,global_seed=1234){
#   message("Running downsample_balanced_iterative")
#   require(dplyr)
#   metadata = as.data.frame(metadata)
#   # group and count full data
#   probs_for_sample = metadata %>% dplyr::ungroup() %>% dplyr::group_by(!!dplyr::sym(predictor_var)) %>% dplyr::add_count(name = "per_group") %>% dplyr::select(Cell_ID,per_group,!!dplyr::sym(predictor_var))
#   probs_for_sample = as.data.frame(probs_for_sample)
#   # first: remove all NAs
#   probs_for_sample = probs_for_sample[!is.na(probs_for_sample[,predictor_var]),]
#   # init
#   sizes = seq(from=nrow(probs_for_sample), to=target_sample_size,by = -stepsize)
#   #print(sizes)
#   n_classes = length(unique(probs_for_sample[,predictor_var]))
#   tmp_samples=probs_for_sample
#   # loop over sizes and remove iteratively
#   cc=0
#   for(size in sizes[1:(length(sizes)-1)]){
#     cc=cc+1
#     # find classes that should be left out
#     which_to_downsample = table(tmp_samples[,predictor_var])>(size/n_classes)
#     # addiotionally add some other classes that are close
#     which_to_downsample = c(which_to_downsample,table(tmp_samples[,predictor_var])>min(table(tmp_samples[,predictor_var])[names(which_to_downsample[which_to_downsample])]))
#     # sample from these two classes
#     set.seed(global_seed)
#     leave_out = sample(tmp_samples$Cell_ID[tmp_samples[,predictor_var] %in% names(which_to_downsample[which_to_downsample])],stepsize)
#     # subset
#     tmp_samples = tmp_samples[!tmp_samples$Cell_ID %in% leave_out,]
#   }
#
#   # return samples IDs
#   return(tmp_samples)
# }
