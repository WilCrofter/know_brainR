
# Return a volume of the BrainWeb phantom which contains its primary visual cortex
# (area 17 or V1.) Note, the returned volume is of type raw.
getArea17 <- function(fname=NA){
  if(is.na(fname) & file.exists("data/subject04_crisp_v.rawb")){
    fname <- "data/subject04_crisp_v.rawb"
  } else if(is.na(fname) & file.exists("../data/subject04_crisp_v.rawb")){
    fname <- "../data/subject04_crisp_v.rawb"
  } else {
    stop("subject04_crisp_v.rawb was not found. Use getArea17(fname) where fname is the path to the file.")
  }
  temp <- readBin(fname, what="raw", n=362*434*362, size=1, signed=FALSE, endian="big")
  dim(temp)<-c(362,434,362)
  temp[(165-19):(165+19), 11:80, (125-19):(125+19)]
}

# The foveal area of the BrainWeb phantom's primary visual cortex appears
# to be ~1 cm (10 voxels) on a side centered around the volume's y axis at x=z=20,
# and consisting of gray matter voxels of minimum x coordinate in this volume.
# This function returns their indices.
fovealRegion <- function(area17){
  ans <- cbind(x_index=rep(seq(20-10, 20+10),20), y_index=NA, z_index=rep(seq(20-10, 20+10),each=20))
  for(n in 1:nrow(ans)){
    ans[n,2] <- min(which(as.integer(area17[ans[n,1],,ans[n,3]]) == 3))
  }
  ans
}

