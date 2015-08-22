
fracAbs <- function(e){
  nx = dim(e$state)[2]
  ny = dim(e$state)[3]
  nz = dim(e$state)[4]
  # Fraction of voxel energies absorbed per step
  frac_abs <- sapply(as.vector(e$state[1,,,]), function(i)e$pAbsorbtion(i))
  dim(frac_abs) <- c(nx,ny,nz)
  frac_abs
}

# Stub for outflow in +x direction.
# TODO: figure out how to control dimension and direction
fracBdry <- function(e){
  nx = dim(e$state)[2]
  ny = dim(e$state)[3]
  nz = dim(e$state)[4]
  # Fraction of voxel photons encountering boundary per step
  ans <- sapply(as.vector(e$state[1,,,]), function(i)e$pBoundary(i))/6
  dim(ans) <- c(nx,ny,nz)
  v1 <- as.vector(e$state[1,-nx,,])
  v2 <- as.vector(e$state[1,-1,,])
  ans[-nx,,] <- ans[-nx,,]*array(mapply(e$pFlow, v1, v2), c(nx-1,ny,nz))
  ans
}

testme <- function(e){
temp <- fracBdry(e)
e$state[1,20:21,40,20]
print(temp[20,40,20])
# [1] 0.02311349
print(e$pBoundary(22)*e$pFlow(22,2))
# [1] 0.1386809
}


voxSim <- function(e, nsteps){
  nx = dim(e$state)[2]
  ny = dim(e$state)[3]
  nz = dim(e$state)[4]
  while(e$step < e$step + nsteps){
    # copying the current state of energy seems necessary
    energy <- e$state[2,,,]
    for(iy in 1:ny){
      for(ix in 1:nx){
        for(iz in 1:nz){
          id <- e$state[1,ix,iy,iz]
          if(id == 0)next # don't care about background
          erg <- energy[ix,iy,iz]
          # Absorption
          ab <- erg*e$pAbsorbtion(id)
          e$state[2,ix,iy,iz] <- e$state[2,ix,iy,iz] - ab
          e$state[3,ix,iy,iz] <- e$state[3,ix,iy,iz] + ab
          # Boundary encounters
          bd <- erg*e$pBoundary(id)/6 # (per boundary)
          # Flow across each boundary TODO: fix
          dx <- c(ix-1, ix+1); dx <- dx[dx >= 1 & dx <= nx]
          dy <- c(iy-1, iy+1); dy <- dy[dy >= 1 & dy <= ny]
          dz <- c(iz-1, iz+1); dz <- dz[dz >= 1 & dz <= nz]
          for(i in dx){
            for(j in dy){
              for(k in dz){
                outflow <- erg*bd*e$pFlow(id, e$state[1,i,j,k])
                e$state[2,ix,iy,iz] <- e$state[2,ix,iy,iz] - outflow
                e$state[2,i,j,k] <- e$state[2,i,j,k] + outflow
              }
            }
          }
          # Treat volume boundaries as absorbing barriers
          if(ix == 1 | ix == nx)e$state[3,ix,iy,iz] <- e$state[3,ix,iy,iz] + bd
          if(iy == 1 | iy == ny)e$state[3,ix,iy,iz] <- e$state[3,ix,iy,iz] + bd
          if(iz == 1 | iz == nz)e$state[3,ix,iy,iz] <- e$state[3,ix,iy,iz] + bd
        }
      }
    }
    # TODO: between steps
    e$step <- e$step + 1
  }
  e
}