
#Given p and u, both 3 dim
#find scalar alpha
findAlph <- function(p,u){
  if (u>=0) alpha <- (1-p)/u
  else alpha <- -p/u
  alpha
}
distance <- function(p,u){
  short <- min(sapply(1:3,function(n){findAlph(p[n],u[n])}))
  short
}
getAll <- function(reps,ranseed){
  set.seed(ranseed)
  p <- matrix(runif(3*reps),reps,3)
  u <- rusphere(reps)
  sapply(1:reps,function(n){distance(p[n,],u[n,])})
  
}
distn <- function(mua,mus,reps,ranseed){
  dist <- getAll(reps,ranseed)
  absorp <- rexp(reps,mua)
  scatter <- rexp(reps,mus)
  pabsorp <- mean(absorp < pmin(scatter,dist))
  pboundary <- mean(dist<pmin(scatter,absorp))
  se_absorp <- sqrt(pabsorp*(1-pabsorp)/reps)
  se_boundary <- sqrt(pboundary*(1-pboundary)/reps) 
  list(pabsorp=pabsorp,pboundary=pboundary,se_absorp=se_absorp,se_boundary=se_boundary)
}
#pr scatter and nothing else
pscatter <- function(mua,mus,reps,ranseed){
  dist <- getAll(reps,ranseed) 
  absorp <- rexp(reps,mua)
  scatter <- rexp(reps,mus)
  prscatt <- mean(scatter<pmin(absorp,dist))
  se_scatt <- sqrt(prscatt*(1-prscatt)/reps)
  list(prscatt=prscatt,se_scatt=se_scatt)
}
pother <- function(mua,mus,reps,ranseed){
  dist <- getAll(reps,ranseed) 
  absorp <- rexp(reps,mua)  
  pabsorb <- mean(absorp<dist)
  pbound <- 1- pabsorb
  se <- sqrt(pabsorb*pbound/reps)
  list(cpabsorb=pabsorb,cpbound=pbound,se=se)
}