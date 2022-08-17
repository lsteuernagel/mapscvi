##########
### reconstruct_levels
##########

#' Reconstruct high cluster levels based on provided edgelist
#'
#' After prediction or label propagation using a low clustering level this function can be used to reconstruct higher clustering levels.
#' Is a wrapper around 'recursive_reconstruction'
#' Requires a clustertree defined in a simple edgelist. User has to ensure that provided input is a valid edgelist!
#'
#' @param edgelist dataframe or matrix with two columns (edgelist: to and from). optional third column giving the 'level' of 'from' which will be used to determine reconstructed levels and name columns
#' @param input_annotation query project name. defaults to 'query'
#' @param level_prefix regex used to extract 'level' if third column is missing
#' @param result_prefix prefix to add to result column names
#'
#' @return dataframe with columns for input_annotation and all clustering levels (above) in edgelist
#'
#' @export
#'
#' @import stringr

reconstruct_levels = function(edgelist,input_annotation,level_prefix="K[0-9]+",result_prefix = ""){

  # TODO: I am not checking if edgelist is truly a valid edgelist!

  # check columns of edgelist and prepare level
  if(ncol(edgelist)==1){
    stop("Error: Please provide a dataframe or matrix with the first two columns correpsonding to 'to' and 'from'")
  }else if(ncol(edgelist)==2){
    colnames(edgelist)[1:2] = c("from","to")
    # infer level
    message("Inferring level with prefix: '",level_prefix,"'")
    edgelist$level = stringr::str_extract(string = edgelist$from,pattern = level_prefix)
    if(all(is.na(edgelist$level))){stop("Error: Cannot infer level. Stopping execution.")}
  }else if(ncol(edgelist)>2){
    # subset to first three and rename third
    edgelist = edgelist[,1:3]
    colnames(edgelist)[1:3] = c("from","to","level")
  }

  # check that input_annotation are part of to
  if(length(setdiff(unique(input_annotation),unique(edgelist$to)))>0){
    stop("Error: Some input annotations are not part of provided edgelist: ",paste0(setdiff(unique(input_annotation),unique(edgelist$to)),collapse = " | "))
  }
  # reconstruct levels
  res= recursive_reconstruction(edgelist,input_annotation)
  colnames(res) = paste0(result_prefix,colnames(res))
  #return
  return(res)
}

##########
### recursive_reconstruction
##########

#' Reconstruct high cluster levels based on provided edgelist
#'
#' Recursively go through edgelist and reconstruct higher-level annotation of input_annotation
#'
#' @param edgelist dataframe or matrix with two columns (edgelist: to and from). optional third column giving the 'level' of 'from' which will be used to determine reconstructed levels and name columns
#' @param input_annotation a vector of labels (per cell). higher level annotations will constrcuted from this
#'
#' @return matrix with columns for input_annotation and all clustering levels (above) in edgelist
#'
#' @export

recursive_reconstruction = function(edgelist,input_annotation){
  current_level = edgelist$level[edgelist$to==input_annotation[1]]
  current_annotation = as.character(sapply(input_annotation,function(x,edgelist){edgelist$from[edgelist$to == x]},edgelist=edgelist))
  if(current_level %in% edgelist$level & !is.na(current_level)){
    next_annotation = recursive_reconstruction(edgelist,current_annotation)
    df = cbind(current_annotation,next_annotation)
    colnames(df)[1]=current_level
  }else{
    df = NULL
  }
  return(df)
}


##########
### add_paired_annotation
##########

#' Return annotation B based on annotation A
#'
#' ..
#'
#' @param input_annotation vector of input annotation A
#' @param reference_mapping dataframe with mapping of annotations A and B. if NULL will be inferred from reference_annotations
#' @param reference_annotations dataframe with all annotations A and B. Will be used to build reference_mapping. Defaults to NULL. Ignored if reference_mapping is provided or NULL.
#'
#' @return vector of same length as input_annotation with annotation B (based on annotation A)
#'
#' @export
#'

add_paired_annotation = function(input_annotation,reference_mapping=NULL,reference_annotations=NULL){

  # if no reference mapping: infer
  if(is.null(reference_mapping)){
    if(!is.null(reference_annotations)){
      reference_mapping = reference_annotations[!duplicated(reference_annotations[,1:2]),]
      rownames(reference_mapping) = reference_mapping[,1]
    }else{
      stop("Error: Please provide either reference_mapping or reference_annotations.")
    }
  }
  # check mapping
  if(max(table(reference_mapping[,1])) > 1 | max(table(reference_mapping[,2])) > 1){
    stop("Error: Please provide a reference_mapping with 1:1 mappings.")
  }
  # check that input_annotation is in reference_mapping
  if(length(setdiff(unique(input_annotation),unique(reference_mapping[,1])))>0){
    stop("Error: Some input annotations are not part of first column of provided reference_mapping: ",paste0(setdiff(unique(input_annotation),unique(reference_mapping[,1])),collapse = " | "))
  }
  # infer annotation B
  new_annotation = as.character(sapply(input_annotation,function(x,reference_mapping){reference_mapping[reference_mapping[,1] == x,2]},reference_mapping=reference_mapping))
  #return vector
  return(new_annotation)
}

