# #####
# ## load neuron map
# #####
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
# reference_hypoMap@assays$RNA@counts <- empty_matrix
# reference_hypoMap@assays$RNA@data <- empty_matrix # error is okay
# # overwrite other slots with empty dummy matrix
# dummy=matrix(data = as.numeric())
# reference_hypoMap@assays$RNA@scale.data <- dummy[,-1] # error is okay
# reference_hypoMap@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
# print(object.size(reference_hypoMap) / 1000000)
# reference_hypoMap_neurons = reference_hypoMap
#
# ## downsample:
# ## source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonize/utils.R") or See function below!!!!
# cluster_downsample = downsample_balanced_iterative(reference_hypoMap_neurons@meta.data,target_sample_size = 20000,predictor_var = "K169_pruned",stepsize = 500,global_seed = 12345)
# table(cluster_downsample$K169_pruned)
#
# reference_hypoMap_neurons = subset(reference_hypoMap_neurons,cells = cluster_downsample$Cell_ID)
# print(object.size(reference_hypoMap_neurons) / 1000000)
#
# save(reference_hypoMap_neurons,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_neurons.RData")
#
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
# reference_hypoMap_full@assays$RNA@counts <- empty_matrix
# reference_hypoMap_full@assays$RNA@data <- empty_matrix # error is okay
# # overwrite other slots with empty dummy matrix
# dummy=matrix(data = as.numeric())
# reference_hypoMap_full@assays$RNA@scale.data <- dummy[,-1] # error is okay
# reference_hypoMap_full@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
# print(object.size(reference_hypoMap_full) / 1000000)
#
# ## downsample:
# ## source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonize/utils.R") or See function below!!!!
# cluster_downsample = downsample_balanced_iterative(reference_hypoMap_full@meta.data,target_sample_size = 20000,predictor_var = "Curated_CellType",
#                                                    stepsize = 500,global_seed = 12345)
# table(cluster_downsample$Curated_CellType)
#
# reference_hypoMap_full = subset(reference_hypoMap_full,cells = cluster_downsample$Cell_ID)
# print(object.size(reference_hypoMap_full) / 1000000)
#
# # save
# save(reference_hypoMap_full,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_full.RData")
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
# #suffix ="mapped_data_yeo_romanov" # a name
# query_romanov_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/mapped_data_Romanov_neurons/mapped_data_Romanov_neurons.h5Seurat" # seurat object to load
# query_romanov = SeuratDisk::LoadH5Seurat(query_romanov_path)
# query_romanov@reductions = list()
# query_romanov@meta.data = query_romanov@meta.data[,1:20]
# dummy=matrix(data = as.numeric())
# query_romanov@assays[["RNA"]]@data <- dummy[,-1] # error is okay
# print(object.size(query_romanov) / 1000000)
# #save testdata
# ##save(query_romanov,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/query_romanov.RData")
# #
# require(scRNAseq)
# sce_lamanno_da <- LaMannoBrainData(which = "mouse-adult",ensembl=FALSE)
# object.size(sce_lamanno_da) / 1000000
# ##save(sce_lamanno_da,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/sce_lamanno_da.RData")
#
# #save testdata
# save(query_romanov,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/query_romanov.RData")
# save(sce_lamanno_da,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/sce_lamanno_da.RData")
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
