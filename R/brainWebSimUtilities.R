#' Convert an array from type raw to type numeric
#' @param raw_array an array of type raw
#' @return an equivalent array of type numeric
raw2numeric <- function(raw_array){
  d <- dim(raw_array)
  raw_array <- as.numeric(raw_array)
  dim(raw_array) <- d
  raw_array
}

#' Return an array of functions to access statistical tables for voxel-level simulation 
#' on BrainWeb volumes.
brainWebAccessors <- function(vox_prob_files, boundary_crossing_files){
  # Load tables into this environment
  readandmerge <- function(files){
    ltabs <- lapply(files, function(f)read.table(f, header=TRUE, sep=',', as.is=TRUE))
    ans <- ltabs[[1]]
    for(tbl in ltabs[-1])ans <- rbind(ans, tbl)
    ans
  }
  vox_probs <- readandmerge(vox_prob_files)
  # Tag vox_probs rows by id for easy lookup
  rownames(vox_probs) <- vox_probs$id
  bdry_probs <- readandmerge(boundary_crossing_files)
  # Remove readandmerge so it doesn't end up in the environment 
  # of the returned functions.
  rm(readandmerge)
  # Function to return absorption probabilities
  pAbsorption <- function(id)if(id==0){1}else{vox_probs[as.character(id),"p_absorb"]}
  # Function to return probability of hitting a boundary
  pBoundary <- function(id)if(id==0){0}else{vox_probs[as.character(id),"p_boundary"]}
  # Function to return probability of flow from id1 to id2
  pFlow <- function(id1, id2){
    if(id1==0)return(0)
    # id's exceeding 11 are gray matter
    if(id1 > 11)id1 <- 2
    if(id2 > 11)id2 <- 2
    bdry_probs[id1, 2 + id2]
  }
  c(pAbsorption=pAbsorption, pBoundary=pBoundary, pFlow=pFlow)
} 

