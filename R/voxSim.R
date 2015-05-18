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