#
# require(devtools)
# require(Seurat)
#
# ## load yeo query
# suffix ="mapped_data_yeo_neurons" # a name
# query_seurat_object_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/mapped_data_yeo_full/mapped_data_yeo_full.h5Seurat" # seurat object to load
# subset_col = "predicted_Curated_Class" # specifiy metdata column to be usedfor subsetting (e.g. subset to neurons when mapping on neuron ref)
# subset_values = "Neurons"
# query_seurat_object = SeuratDisk::LoadH5Seurat(query_seurat_object_path)
#
# ## load romanov query
# suffix ="mapped_data_yeo_romanov" # a name
# query_seurat_object_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/mapped_data_Romanov_neurons/mapped_data_Romanov_neurons.h5Seurat" # seurat object to load
# query_seurat_object = SeuratDisk::LoadH5Seurat(query_seurat_object_path)
# query_seurat_object@reductions = list()
# query_seurat_object@meta.data = query_seurat_object@meta.data[,1:20]
#
# #save testdata
# save(query_seurat_object,
#      file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/query_test_object.RData")
#
# # set model
# map_name = "hypothalamus_neurons_reference" # reference map
# map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/"
# model_path = paste0(map_path,map_name,"_model/")
#
# ## load neuron map
# map_name = "hypothalamus_neurons_reference"
# map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/"
# map_seurat_path = paste0(map_path,map_name,".h5Seurat")
# neuron_map_seurat = SeuratDisk::LoadH5Seurat(map_seurat_path)
#
# #save
# reference_map_reduc = neuron_map_seurat@reductions$scvi
# reference_map_umap = neuron_map_seurat@reductions$umap_scvi
# reference_map_metadata = neuron_map_seurat@meta.data
# save(reference_map_reduc,reference_map_umap,reference_map_metadata,
#                            file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/hypoMap_dimred.RData")
#
#
#
#

# reference_hypoMap = CreateSeuratObject(neuron_map_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = neuron_map_seurat@meta.data)
# reference_hypoMap@reductions = neuron_map_seurat@reductions
# dummy=matrix(data = as.numeric())
# reference_hypoMap@assays$RNA@counts <- dummy[,-1] # error is okay
# reference_hypoMap@assays$RNA@data <- dummy[,-1] # error is okay
# reference_hypoMap@assays$RNA@scale.data <- dummy[,-1] # error is okay
# reference_hypoMap@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
# save(reference_hypoMap,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/tmp_mapscvi/reference_hypoMap.RData")
