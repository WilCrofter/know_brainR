voxSim <- function(e, nsteps){
  while(e$step < e$step + nsteps){
    # copying the current state of energy seems necessary
    energy <- e$state[1,,,]
    nx = dim(e$state)[2]
    ny = dim(e$state)[3]
    nz = dim(e$state)[4]
    for(ix in 1:nx){
      for(iy in 1:ny){
        for(iz in 1:nz){
          id <- e$state[1,ix,iy,iz]
          erg <- energy[ix,iy,iz]
          # Absorption
          ab <- erg*e$pAbsorbtion()
          e$state[2,ix,iy,iz] <- e$state[2,ix,iy,iz] - ab
          e$state[3,ix,iy,iz] <- e$state[3,ix,iy,iz] + ab
          # Boundary encounters
          bd <- erg*e$pBoundary(id)/6 # (per boundary)
          # Flow across each boundary
          dx <- max(ix-1,1):min(ix+1,nx)
          dy <- max(iy-1,1):min(iy+1,ny)
          dz <- max(iz-1,1):min(iz+1,nz)
          for(i in dx){
            for(j in dy){
              for(k in dz){
                outflow <- erg*bd*e$pFlow(id, 3$state[1,i,j,k])
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
    e$step <- e$step + 1
  }
  e
}