#' A second proposal for a voxel-level simulator
#' 
#' It needs a state array for the phantom region of interest,
#' which would be 4D array, the dimension holding the state
#' and the remaining 3 identifying the voxel.  This prototype 
#' uses a state consisting of 3 numbers: tissue, energy, and cum_absorbed.
#' Thus state[1, i, j, k] is the tissue identifier of voxel i,j,k,
#' state[2, i, j, k] is its internal energy and state[3, i, j, k]
#' is its cumulative energy absorbed in chromophores and dyes. 

# This function returns a dummy state array
dummy_state <- function(xdim, ydim, zdim){
  dummy <- rep(0, 3*xdim*ydim*zdim)
  # randomly assigned tissues
  dummy[seq(1, 3*xdim*ydim*zdim, by=3)] <- sample(1:11, xdim*ydim*zdim, replace=TRUE)
  # dimension the array
  dim(dummy) <- c(3, xdim, ydim, zdim)
  # name the rows (first dimension)
  row.names(dummy) <- c("tissue", "energy", "cum_absorbed")
  # just for illustration, make the xy plane at z=2 a source
  dummy[2,,,2] <- 1
  dummy
}

# Given 3D tissue-type arrays of the same shape, M1 and M2, and
# a table of source-to-target boundary crossing probabilities,
# return an array of probabilities having the same shape as M1 and M2.
# TODO: optimize
flow_fractions <- function(M1, M2, bdry_probs){
  ans <- numeric(length(M1))
  d <- dim(ans) <- dim(M1)
  for(i in 1:(d[1])){
    for(j in 1:(d[2])){
      for(k in 1:d[3]){
        ans[i,j,k] <- bdry_probs[M1[i,j,k], 2+M2[i,j,k]]
      }
    }
  }
  ans
}

# An "optimized" form of flow_fractions
# NOTE: lightly tested, but seems to get the same answers as flow_fractions
flow_fractions_opt <- function(M1, M2, bdry_probs){
  ans <- mapply(function(i,j)bdry_probs[i,j+2], as.vector(M1), as.vector(M2))
  dim(ans) <- dim(M1)
  ans
}

#' The following allegedly does 1 step. It surely needs improvement. For instance, we
#' wouldn't make state_array an explicit parameter because it's big and would be copied.
#' We might want more than 1 step. The absorbed and boundaries arrays are big, so we
#' may have to do the associated operations in pieces, e.g., by looping on 1 or 2
#' dimensions.
#' vox_probs and bdry_probs were read in from data/vox_probs.csv and
#' data/boundary_crossing_probs.csv which we created the other day.
#' The function is supposed to conserve energy so
#'     sum(state_array[c("energy", "cum_absorbed"),,,])
#' should be the same before and after a step.
vox_sim <- function(state_array, vox_probs, bdry_probs){
  xdim <- dim(state_array)[2]
  ydim <- dim(state_array)[3]
  zdim <- dim(state_array)[4]
  #   Calculate absorption, add to cumulative absorption
  #   Recall that component 2 of a voxel's state is energy and component 1 is tissue.
  absorbed <- state_array[2,,,]*vox_probs[state_array[1,,,], "p_absorb"]
  state_array["cum_absorbed",,,] <- state_array["cum_absorbed",,,] + absorbed
  #   Calculate total energy hitting EACH boundary (hence the 1/6)
  boundaries <- (1/6)*state_array[2,,,]*vox_probs[state_array[1,,,], "p_boundary"]
  #   Subtract absorbed from internal energy, conserving the total
  state_array[2,,,] <- state_array[2,,,] - absorbed
  rm(absorbed) # free memory
  #   Calculate flow to each neighbor, subtract from source's and add to neighbor's energy
  # TODO: In general we'll NOT want tissue type 0 to be a source. I think we might
  # finesse this by giving it and things like detectors absorption probs of 1, so
  # that the corresponding boundary (and scattering) probs will be 0.
  # 
  # +x direction. Recall we don't want energy flowing from voxels at the edge of the region
  # e.g., from voxels having x=1 or x=xdim. So, e.g., flow in the +x direction would be
  # from voxel 2,j,k  to voxel 3,j,k ... voxel dimx-1, j, k to voxel dimx
  flow <- boundaries[2:(xdim-1),,] * flow_fractions(state_array[1, 2:(xdim-1),,], 
                                                     state_array[1, 3:xdim,,], bdry_probs) 
  #      subtract from sources
  state_array[2, 2:(xdim-1),,] <- state_array[2, 2:(xdim-1),,] - flow
  #      add to +x neighbors, conserving energy
  state_array[2, 3:xdim,,]     <- state_array[2,3:xdim,,]      + flow
  #
  # -x direction: sources are now 2,j,k ... (xdim-1), j, k
  #               sinks   are     1,j,k ... (xdim-2), j, k
  flow <- boundaries[2:(xdim-1),,] * flow_fractions(state_array[1, 2:(xdim-1),,], 
                                                     state_array[1, 1:(xdim-2),,], bdry_probs) 

  #      subtract from sources
  state_array[2, 2:(xdim-1),,] <- state_array[2, 2:(xdim-1),,]    - flow
  #      add to -x neighbor, conserving energy
  state_array[2, 1:(xdim-2),,]     <- state_array[2,1:(xdim-2),,] + flow   
  #
  # + y direction
  flow <- boundaries[,2:(ydim-1),] * flow_fractions(state_array[1,,2:(ydim-1),], 
                                                     state_array[1,,3:ydim,], bdry_probs) 
  #      subtract from sources
  state_array[2,,2:(ydim-1),] <- state_array[2,,2:(ydim-1),] - flow
  #      add to +y neighbors, conserving energy
  state_array[2,,3:ydim,]     <- state_array[2,,3:ydim,]      + flow
  # -y direction
  flow <- boundaries[,2:(ydim-1),] * flow_fractions(state_array[1,,2:(ydim-1),], 
                                                    state_array[1,,1:(ydim-2),], bdry_probs) 
  #      subtract from sources
  state_array[2,,2:(ydim-1),] <- state_array[2,,2:(ydim-1),]    - flow
  #      add to -y neighbor, conserving energy
  state_array[2,,1:(ydim-2),]     <- state_array[2,,1:(ydim-2),] + flow
  #
  # +z direction
  flow <- boundaries[,,2:(zdim-1)] * flow_fractions(state_array[1,,,2:(zdim-1)], 
                                                    state_array[1,,,3:zdim], bdry_probs) 
  #      subtract from sources
  state_array[2,,,2:(zdim-1)] <- state_array[2,,,2:(zdim-1)] - flow
  #      add to +y neighbors, conserving energy
  state_array[2,,,3:zdim]     <- state_array[2,,,3:zdim]      + flow
  # -z direction
  flow <- boundaries[,,2:(zdim-1)] * flow_fractions(state_array[1,,,2:(zdim-1)], 
                                                    state_array[1,,,1:(zdim-2)], bdry_probs) 
  #      subtract from sources
  state_array[2,,,2:(zdim-1)] <- state_array[2,,,2:(zdim-1)]    - flow
  #      add to -y neighbor, conserving energy
  state_array[2,,,1:(zdim-2)]     <- state_array[2,,,1:(zdim-2)] + flow  
  #
  state_array
}


