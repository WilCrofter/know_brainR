## Utility to prepare BrainWeb phantom Area 17 for voxel-level simulation.

#' Return a volume of the BrainWeb phantom which contains its primary visual cortex
#' (area 17) with foveal area properly stained.
#' @param fname_BrainWeb_Phantom path to the BrainWeb file
#' @param fovealTissueID an integer between 21 and 34 (inclusive) indicating level of dye concentration. 
getBrainWebArea17 <- function(fname_BrainWeb_Phantom, fovealTissueID){
  a17 <- readBin(fname_BrainWeb_Phantom, what="raw", 
                  n=362*434*362, size=1, signed=FALSE, endian="big")
  dim(a17)<-c(362,434,362)
  a17 <- a17[(165-19):(165+19), 11:80, (125-19):(125+19)]
  # The foveal area of the BrainWeb phantom's primary visual cortex appears
  # to be ~1 cm (10 voxels) on a side centered around the volume's y axis at x=z=20,
  # and consisting of gray matter voxels of minimum x coordinate in this volume.
  for(ix in seq(20-10, 20+10)){
    for(iz in seq(20-10, 20+10)){
      a17[ix, min(which(as.integer(a17[ix,,iz])==3))] <- fovealTissueID
    }
  }
  a17
}
