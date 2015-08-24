source("R/voxSim.R")
source("R/brainWebSimUtilities.R")
source("R/area17.R")

absorbedStained <- function(e){
  # initialize if necessary
  if(!exists("stained_area", envir=e, inherits=FALSE)){
    e$absorbed_stained <- numeric()
    e$stained_area <- matrix(0,3,0)
    for(ix in seq(20-5, 20+5)){
      for(iz in seq(20-5, 20+5)){
        e$stained_area <- 
          cbind(e$stained_area,
                c(ix, which(e$state[1,ix,,iz] == e$fovealTissueID), iz))
      }
    }
  }
  # compute currrent value
  cval <- sum(apply(e$stained_area, 2, function(v){e$state[3,v[1],v[2],v[3]]}))
  # append current value to record
  e$absorbed_stained <- c(e$absorbed_stained, cval)
}

impulse <- function(e)laserExcitationForArea17(e,1)

btw_fcts <- c(impulse, absorbedStained)

ans <- matrix(0,0,2)

for(ic in 21:34){
  e <- area17env("data", ic, btw_fcts)
  voxSim(e,5000)
  ans <- rbind(ans, c(ic, max(e$absorbed_stained)))
  print(ic)
}

saveRDS(ans, file="ans.rda")