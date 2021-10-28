#####
## load neuron map
#####
map_name = "hypothalamus_neurons_reference"
map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/"
map_seurat_path = paste0(map_path,map_name,".rds")
neuron_map_seurat = readRDS(map_seurat_path)

# # make new seurat object with reduced size
reference_hypoMap = CreateSeuratObject(neuron_map_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = neuron_map_seurat@meta.data)
# add reductions back in
reference_hypoMap@reductions = neuron_map_seurat@reductions
# overwrite counts with empty sparse matrix to save space
empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = reference_hypoMap@assays$RNA@counts@Dim, dimnames = reference_hypoMap@assays$RNA@counts@Dimnames)
empty_matrix <- as(empty_matrix, "dgCMatrix")
reference_hypoMap@assays$RNA@counts <- empty_matrix
reference_hypoMap@assays$RNA@data <- empty_matrix # error is okay
# overwrite other slots with empty dummy matrix
dummy=matrix(data = as.numeric())
reference_hypoMap@assays$RNA@scale.data <- dummy[,-1] # error is okay
reference_hypoMap@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
print(object.size(reference_hypoMap) / 1000000)
reference_map_reduc = neuron_map_seurat@reductions$scvi
object.size(reference_map_reduc) / 1000000
object.size(reference_map_metadata) / 1000000
reference_map_umap = neuron_map_seurat@reductions$umap_scvi
reference_map_metadata = neuron_map_seurat@meta.data
object.size(reference_map_metadata) / 1000000
object.size(reference_map_umap) / 1000000
reference_hypoMap_neurons=reference_hypoMap
save(reference_hypoMap_neurons,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_neurons.RData")

# set model
model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_models/hypothalamus_neurons_reference_model/"
system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))

#####
## load full map
#####

project_name = "hypothalamus_full_map"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/"
seurat_file_name = paste0(project_path,project_name,".h5Seurat")
full_map_seurat = SeuratDisk::LoadH5Seurat(seurat_file_name)

# # make new seurat object with reduced size
reference_hypoMap_full = CreateSeuratObject(full_map_seurat@assays$RNA@counts,project = "HypoMap_reference",meta.data = full_map_seurat@meta.data)
# add reductions back in
reference_hypoMap_full@reductions = full_map_seurat@reductions
# overwrite counts with empty sparse matrix to save space
empty_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = reference_hypoMap_full@assays$RNA@counts@Dim, dimnames = reference_hypoMap_full@assays$RNA@counts@Dimnames)
empty_matrix <- as(empty_matrix, "dgCMatrix")
reference_hypoMap_full@assays$RNA@counts <- empty_matrix
reference_hypoMap_full@assays$RNA@data <- empty_matrix # error is okay
# overwrite other slots with empty dummy matrix
dummy=matrix(data = as.numeric())
reference_hypoMap_full@assays$RNA@scale.data <- dummy[,-1] # error is okay
reference_hypoMap_full@assays$RNA@meta.features <- as.data.frame(dummy[,-1])
print(object.size(reference_hypoMap_full) / 1000000)
save(reference_hypoMap_full,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/reference_hypoMap_full.RData")

# set model
model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_models/scVI_hypothalamus_full_map_model/"
system(paste0("cp -r ",model_path," /beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/inst/extdata/models/"))

#####
## load example data
#####

## load romanov query
#suffix ="mapped_data_yeo_romanov" # a name
query_romanov_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/mapped_data_Romanov_neurons/mapped_data_Romanov_neurons.h5Seurat" # seurat object to load
query_romanov = SeuratDisk::LoadH5Seurat(query_romanov_path)
query_romanov@reductions = list()
query_romanov@meta.data = query_romanov@meta.data[,1:20]
#save testdata
save(query_romanov,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/query_romanov.RData")
#
require(scRNAseq)
sce_lamanno_da <- LaMannoBrainData(which = "mouse-adult",ensembl=FALSE)
object.size(sce_lamanno_da)
object.size(sce_lamanno_da) / 1000000
save(sce_lamanno_da,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/sce_lamanno_da.RData")
dummy=matrix(data = as.numeric())
query_romanov@assays[["RNA"]]@data <- dummy[,-1] # error is okay
#save testdata
save(query_romanov,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/query_romanov_object.RData")
save(sce_lamanno_da,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/mapscvi/data/sce_lamanno_da.RData")
