#' Convert an array from type raw to type numeric
#' @param raw_array an array of type raw
#' @return an equivalent array of type numeric
raw2numeric <- function(raw_array){
  d <- dim(raw_array)
  raw_array <- as.numeric(raw_array)
  dim(raw_array) <- d
  raw_array
}

#' Prepare statistical tables for voxel-level simulation on BrainWeb volumes
brainWebTables <- function(vox_prob_files, boundary_crossing_files){
  readandmerge <- function(files){
    ltabs <- lapply(files, function(f)read.table(f, header=TRUE, sep=',', as.is=TRUE))
    ans <- ltabs[[1]]
    for(tbl in ltabs[-1])ans <- rbind(ans, tbl)
    ans
  }
  vox_probs <- readandmerge(vox_prob_files)
  bdry_probs <- readandmerge(boundary_crossing_files)
  # Boundary probs for stained gray matter (id > 20) will be the same as for
  # gray matter itself. Augment tables to reflect this
  # First determine id's for stained gray matter
  idx <- vox_probs$id > 20
  ids <- vox_probs$id[idx]
  tissues <- vox_probs$tissue[idx]
  # append repetitions of gray matter columns for each stain identifier
  for(n in 1:length(ids))bdry_probs[,tissues[n]] <- bdry_probs[, "Gray.Matter"]
  # append repetitions of gray matter row for each stain identifier
  gmr <- bdry_probs[bdry_probs$source_tissue == "Gray Matter",]
  for(n in 1:length(ids)){
    gmr[1,1] <- tissues[n]
    bdry_probs <- rbind(bdry_probs, gmr)
  }
  list(vox_probs=vox_probs, bdry_probs=bdry_probs)
} 

