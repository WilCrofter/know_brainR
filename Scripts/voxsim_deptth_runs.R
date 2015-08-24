# Tracks energy density and absorption as a function of depth in a column beneath
# the area of excitation at the scalp. This file merely sets up the run. Iterate
# for 5000 steps after sourcing.

# Source this file from the parent directory of these files
source("R/voxSim.R")
source("R/brainWebSimUtilities.R")
source("R/area17.R")

#  a function to illuminate an area of the scalp on the first step only
impulse <- function(e)laserExcitationForArea17(e,1)

#  a function to calculate absorption and energy density as a function of depth 
# in a column beneath the excitation
depth <- function(e, how_often=100){
  if(e$step==0 | e$step %% how_often != 0)return()
  if(!exists("erg", envir = e, inherits = FALSE))e$erg <- numeric()
  if(!exists("absr", envir = e, inherits = FALSE))e$absr <- numeric()
  side <- (20-5):(20+5)
  ny <- length(e$state[1,1,,1])
  erg <- absr <- numeric()
  for(iy in 7:ny){
    erg <- c(erg, sum(e$state[2,side,iy,side]))
    absr <- c(absr, sum(e$state[3,side,iy,side]))
  }
  e$erg <- cbind(e$erg, erg)
  e$absr <- cbind(e$absr, absr)
}

phantom <- area17env("data", 28, c(impulse, function(e)depth(e,500)))

voxSim(phantom, 1)
