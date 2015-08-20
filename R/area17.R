## Utilities for voxel-level simulation on BrainWeb phantom Area 17.

source("R/brainWebSimUtilities.R", local=TRUE)

#' Returns an initialized state array for a volume of the BrainWeb phantom 
#' which includes its primary visual cortex (area 17) with foveal area properly stained.
#' @param fname_BrainWeb_Phantom path to the BrainWeb file
#' @param fovealTissueID an integer between 21 and 34 (inclusive) indicating level of dye concentration. 
brainWebArea17 <- function(fname_BrainWeb_Phantom, fovealTissueID){
  a17 <- readBin(fname_BrainWeb_Phantom, what="raw", 
                  n=362*434*362, size=1, signed=FALSE, endian="big")
  dim(a17)<-c(362,434,362)
  a17 <- a17[(165-19):(165+19), 11:80, (125-19):(125+19)]
  # The foveal area of the BrainWeb phantom's primary visual cortex appears
  # to be ~1 cm (10 voxels) on a side centered around the volume's y axis at x=z=20,
  # and consisting of gray matter voxels of minimum x coordinate in this volume.
  fovealTissueID <- as.raw(fovealTissueID)
  for(ix in seq(20-5, 20+5)){
    for(iz in seq(20-5, 20+5)){
      a17[ix, min(which(as.integer(a17[ix,,iz])==3)), iz] <- fovealTissueID
    }
  }
  d <- dim(a17)
  ans <- numeric(3*prod(dim(a17)))
  dim(ans) <- c(3, dim(a17))
  ans[1,,,] <- raw2numeric(a17)
  ans
}

#' Return the indices of the voxels on the scalp surface within
#' an 11x11 square area surrounding voxel x=z=20 of the area17 volume.
#' @param e an environment such that e$state is the state of area 17
#' @return a 121x3 array of scalp coordinates
scalpArea17 <- function(e){
  ans <- integer()
  for(ix in seq(20-5, 20+5)){
    for(iz in seq(20-5, 20+5)){
      iy <- 1+max(which(e$state[1,ix,,iz]==0))
      ans <- rbind(ans, c(ix, iy, iz))
    }
  }
  ans
}

#' A utility to provide excitation at the scalp for simulations involving
#' area 17. This utility itself should not be provided, since an excitation
#' function must have only one parameter--the environment containing the state.
#' Instead, provide a function such as function(e)excitationForArea17(e,500),
#' i.e., a derivative function where the parameter n is fixed.
#' @param e environment containing the state
#' @param n the number of steps for which unit excitation should be provided.
#' @return none. The function has a side effect upon e$state.
excitationForArea17 <- function(e, n){
  if(!exists(scalp_indices, envir=e, inherits=FALSE)){
    e$scalp_indices <- scalpArea17(e)
  }
  if(e$step <= n){
    for(i in 1:nrow(e$scalp_indices)){
      e$state[2, scalp_indices[i,1], scalp_indices[i,2], scalp_incices[i,3]] <-
        1 + e$state[2, scalp_indices[i,1], scalp_indices[i,2], scalp_incices[i,3]]
    }
  }
  invisible()
}
  