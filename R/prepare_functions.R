
# This file contains functions to prepare a query seurat for mapping with mapscvi.

##########
### transform_seurat
##########

#' Update gene ids
#'
#' Helper function to prepare seurat object gene ids for format used in reference
#' Currently does nothing and just returns the same object.
#'
#' @param query_seurat_object Seurat object
#' @param suffix project name to clearly label various steps. defaults to 'query'
#' @param global_seed seed
#'
#' @return formatted seurat object
#'
#' @export
#'
#' @import SeuratObject Seurat
#'
#' @examples

transform_seurat = function(query_seurat_object,suffix="query",global_seed=12345){

  #TODO: check that Cell_ID exists and combine with suffix

  #TODO: procedures to map to gene names used in reference query

  # from human ens ids

  # from human gene names

  # from mouse ens ids

  return(query_seurat_object)

}


##########
### prepare_query
##########

#' Prepare seurat object for mapping functions
#'
#' Generic Function to prepare a Seurat object for embedding
#'
#' @param query_seurat_object Seurat object
#' @param suffix project name to clearly label various steps. defaults to 'query'
#' @param assay which assay from query_seurat_object. defaults to RNA
#' @param subset_col character value: metadata column in query_seurat_object that will be used to subset data. defaults to ''
#' @param subset_values character vector: which values from subset_col should be kept. defaults to NULL which will prevent any subsetting
#' @param normalize boolean: run standard normalization on assay. defaults to TRUE
#' @param global_seed seed
#'
#' @return updated seurat object
#'
#' @export
#'
#' @import Seurat SeuratObject rlang dplyr
#'
#' @examples

prepare_query = function(query_seurat_object,suffix="query",assay="RNA",subset_col="",subset_values=NULL,normalize=TRUE,global_seed=12345){

  # subset if wanted
  if(subset_col %in% colnames(query_seurat_object@meta.data) & !is.null(subset_values)){
    bef = ncol(query_seurat_object)
    query_seurat_object = subset(query_seurat_object,subset = !!rlang::sym(subset_col) %in% subset_values)
    message("Subsetting query object to ",ncol(query_seurat_object)," cells from ",bef," cells" )
  }
  #normalize
  if(dim(neuron_map_seurat@assays$RNA@data)[2] != dim(neuron_map_seurat@assays$RNA@counts)[2]){normalize=TRUE}
  if(normalize){
    message("Normalizing data")
    query_seurat_object <- Seurat::NormalizeData(object = query_seurat_object,assay = assay, verbose = F,normalization.method = "LogNormalize",scale.factor = 10000)
  }

  ## clean up
  query_seurat_object@project.name = suffix
  query_seurat_object@reductions = list()
  #query_seurat_object@misc = list()
  dummy=matrix(data = as.numeric())
  query_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

  return(query_seurat_object)
}

##########
### prepare_query_hypoMap
##########

#' Prepare seurat object for mapping onto hypoMap
#'
#' Function to prepare a Seurat object for embedding with releavnt covariates and metadata columns to work with HypoMap
#'
#' @inheritParams prepare_query
#' @param covariates named vector: names of covariates for hypoMap scvi model that will be added via this function. Changing this parameter or supplying the variables without the proper format can cause probelms! defaults to:
#' @param sex_var column name with annotation of sample sex
#'
#' @return updated seurat object
#'
#' @export
#'
#' @import Seurat SeuratObject rlang dplyr
#'
#' @examples

prepare_query_hypoMap = function(query_seurat_object,suffix="query",assay="RNA",subset_col="",subset_values=NULL,normalize=TRUE,sex_var = "Sex",
                                 covariates=c(batch_var = "Batch_ID",inferred_sex = "inferred_sex.x",rpl_signature_expr_median = "rpl_signature_expr_median"),global_seed=12345){
  # subset if wanted
  if(subset_col %in% colnames(query_seurat_object@meta.data) & !is.null(subset_values)){
    bef = ncol(query_seurat_object)
    query_seurat_object = subset(query_seurat_object,subset = !!rlang::sym(subset_col) %in% subset_values)
    message("Subsetting query object to ",ncol(query_seurat_object)," cells from ",bef," cells" )
  }
  #normalize
  if(dim(neuron_map_seurat@assays$RNA@data)[2] != dim(neuron_map_seurat@assays$RNA@counts)[2]){normalize=TRUE}
  if(normalize){
    message("Normalizing data")
    query_seurat_object <- Seurat::NormalizeData(object = query_seurat_object,assay = assay, verbose = F,normalization.method = "LogNormalize",scale.factor = 10000)
  }

  # add missing variables for hypomap
  batch_var = covariates["batch_var"]
  inferred_sex = covariates["inferred_sex"]
  rpl_signature_expr_median = covariates["rpl_signature_expr_median"]
  ## check if batch var exists
  if(!batch_var %in% colnames(query_seurat_object@meta.data)){
    message("Batch variable '",batch_var,"' not found. Adding with value '",paste0(suffix,"_batch_1"),"'")
    query_seurat_object@meta.data[,batch_var] = paste0(suffix,"_batch_1")
  }

  ## check if inferred sex exists
  if(!inferred_sex %in% colnames(query_seurat_object@meta.data)){
    message("Variable '",inferred_sex,"' not found. Adding to data","'")
    xist_expr <- tryCatch({
      xist_expr = Seurat::FetchData(query_seurat_object,vars = "Xist")$Xist
      xist_expr[xist_expr>0] = 1
      xist_expr
    },
    error=function(cond) {
      message("Cannot determine Xist expression level")
      # Choose a return value in case of error
      return(NA)
    })
    if(is.na(xist_expr[1])){
      if(sex_var %in% colnames(query_seurat_object@meta.data)){
        if(c("M") %in% unique(query_seurat_object@meta.data[,sex_var])){
          message("Cannot use Xist. Using '",sex_var,"' to infer")
          query_seurat_object@meta.data$inferred_sex.x = query_seurat_object@meta.data[,sex_var]
        }
      }else{
        message("Impossible to infer Sex. Setting all cells to 'M'")
        query_seurat_object@meta.data$inferred_sex.x = "M"
      }
    }else{
      if(sex_var %in% colnames(query_seurat_object@meta.data)){
        if(c("M") %in% unique(query_seurat_object@meta.data[,sex_var])){
          message("Found '",sex_var,"' annotation. Using custom procedure to compare with Xist expression")
          tmp_df = data.frame(Cell_ID = query_seurat_object@meta.data$Cell_ID, Sample_ID = query_seurat_object@meta.data$Sample_ID,xist_expr,Sex=query_seurat_object@meta.data[,sex_var]) %>% group_by(Sample_ID) %>%
            dplyr::add_count(name="n_cells") %>% dplyr::mutate(n_xist = sum(xist_expr)) %>% dplyr::mutate(xist_pct = n_xist / n_cells)
          tmp_df$inferred_sex.x = tmp_df$xist_expr
          tmp_df$inferred_sex.x[tmp_df$Sex=="F" & tmp_df$xist_pct > 0.5] = "1" # overwrite clear samples (even if there are some Xist outlier cells
          tmp_df$inferred_sex.x[tmp_df$Sex=="M" & tmp_df$xist_pct < 0.2] = "0" # overwrite clear samples (even if there are some Xist outlier cells)
          tmp_df$inferred_sex.x[tmp_df$inferred_sex.x=="0"] = "M"
          tmp_df$inferred_sex.x[tmp_df$inferred_sex.x=="1"] = "F"
          tmp_meta = left_join(query_seurat_object@meta.data,tmp_df[,c("Cell_ID","inferred_sex.x")],by="Cell_ID")
          rownames(tmp_meta) = tmp_meta$Cell_ID
          query_seurat_object@meta.data = tmp_meta
        }else{
          message("Setting inferred sex based on Xist expression")
          query_seurat_object@meta.data$inferred_sex.x[xist_expr>0] = "F"
          query_seurat_object@meta.data$inferred_sex.x[is.na(query_seurat_object@meta.data$inferred_sex.x)] = "M"
        }
      }else{
        message("Setting inferred sex based on Xist expression")
        query_seurat_object@meta.data$inferred_sex.x[xist_expr>0.1] = "F"
        query_seurat_object@meta.data$inferred_sex.x[is.na(query_seurat_object@meta.data$inferred_sex.x)] = "M"
      }
    }
  }
  ## check if rpl sig exists
  if(!rpl_signature_expr_median %in% colnames(query_seurat_object@meta.data)){
    message("Variable '",rpl_signature_expr_median,"' not found. Adding to data.")
    rpl_signature = as.character(c("Rpl32", "Rpl26", "Rpl22l1", "Rps19", "Rpl39", "Rps9", "Rps15", "Rpl27a", "Rps25", "Rps20", "Rpl36al", "Rps8", "Rpl6", "Rps23", "Rpl19", "Cox7a2l", "Rpl7", "Rps26-ps1", "Rps27a", "Rps16", "Rpl29", "Rps3a1", "Rps4x", "Rps6", "Rpl17", "Rps18-ps3", "Rpl9-ps6", "Rps24", "Rps12-ps3", "Rps7", "Rpl21", "Rpl36-ps3", "Rpl13-ps3", "Rpl23a-ps3", "Rpl10a", "Rpl6l", "Rpl27-ps3", "Rpl7a-ps5", "Cox7c", "Jund", "Rpl10"))
    rpl_signature_expr = Seurat::FetchData(query_seurat_object,vars = rpl_signature)
    query_seurat_object@meta.data$rpl_signature_expr_median = apply(rpl_signature_expr,1,median,na.rm=TRUE)
  }


  ## clean up
  query_seurat_object@project.name = suffix
  query_seurat_object@reductions = list()
  #query_seurat_object@misc = list()
  dummy=matrix(data = as.numeric())
  query_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

  return(query_seurat_object)
}

