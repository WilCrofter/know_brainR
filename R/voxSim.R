#' A second proposed outline for a voxel-level simulator
#' 
#' First, prepare a state array for the phantom region of interest:
#'    4D array: i,j,k to id voxel and 4-long state: 
#'    tissue, #photons, #incoming, cumulative absorption
#'
vox_sim <- function(n, state_array, vox_probs, bdry_probs){
  # First pass:
  #   Calculate absorption, subtract from #photons, add to cumulative absorption
  #   Calculate flow to each neighbor, subtract from self, add to neighbor's incoming
  # Second pass:
  #   Add #incoming to #photons and reset #incoming to 0.
}

# Fake phantom
phantom <- rep(0, 4^4)
phantom[seq(1, 4^4, by=4)] <- sample(1:11, 4^3, replace=TRUE)
dim(phantom) <- c(4,4,4,4)
row.names(phantom) <- c("tissue", "n", "nin", "abs")
phantom
# [,1] [,2] [,3] [,4]
# tissue    5   11    8    2
# n         0    0    0    0
# nin       0    0    0    0
# abs       0    0    0    0
# ...
# Thus phantom[, i, j, k] is the state of voxel i, j, k
phantom[,3, 2, 4]
# tissue      n    nin    abs 
#      2      0      0      0

#' A proposed outline for a voxel-level simulator. Needs sanity checks.
#' 
#' @param n number of steps to take
#' @param voxel_states a list giving coordinates, i, j, k, internal energy, and other voxel state info.
#' @param excitation a function of n and voxel_states which adds energy to certain voxels in voxel_states list
#' @param get_voxel_type a function of voxel index i, j, k which may return special values such as boundary
vox_sim <- function(n, voxel_states, excitation, get_voxel_types){
  # For 1:n
  # First pass:
  # For each voxel in voxel_states
  #   1. Calculate the proportion of its energy lost to absorption and subtract it from internal energy.
  #   2. Add the absorbed proportion to a running total.
  #   3. Calculate the proportion of energy lost across each of its 6 boundaries and subtract
  #      the total from its internal energy. Save the 6 values for 2nd pass.
  # Second pass:
  # For each voxel in voxel states, add the energy lost across boundaries to each of its 6
  #   neighbors, adding the neighbor to voxel_states if necessary, unless it is a boundary
  #   voxel.
  # end for loop
  # return voxel_states
}
