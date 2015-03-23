#computes probability of reflection
#given angle of incidence (in radians)
#and two indices of refraction, 
#source and dest media respectively
prob_refl <- function(theta1,n1,n2){
  prob <- rep(1,length(theta1))
  cos2 <- rep(0,length(theta1))
  sin1 <- sin(theta1)
  cos1 <- sqrt(1-sin1*sin1)
  sin2 <- (n1/n2) * sin1
  idx <- sin2<=1
  cos2[idx] <- sqrt(1-sin2[idx]^2)
  t1 <-  (abs((n1[idx]*cos1[idx]-n2[idx]*cos2[idx])/(n1[idx]*cos1[idx]+n2[idx]*cos2[idx])))^2
  t2 <-  (abs((n1[idx]*cos2[idx]-n2[idx]*cos1[idx])/(n1[idx]*cos2[idx]+n2[idx]*cos1[idx])))^2  
  prob[idx] <- pmin(1, .5 * (t1+t2))
  prob
}
#compute probability of reflection for angles moving from
#gray matter to white matter
Rgraytowhite <- function(theta){
  n_gray <- dry - (dry-wet)*rep(water_gray,length(theta))
  n_white <- dry - (dry-wet)*rep(water_white,length(theta) )
  prob_refl(theta,n_gray,n_white)
}
#compute angles of refraction in white matter
#given incidence angles coming in from gray matter
refraction_in_white <- function(theta){
  reflprobs <- Rgraytowhite(theta)
  idx_refl <- reflprobs==1.0
  not_reflected <- theta[!idx_refl]
  n_gray <- dry - (dry-wet)*water_gray
  n_white <- dry - (dry-wet)*water_white
  sin2 <- (n_gray/n_white) * sin(not_reflected)
  cos2 <- sqrt(1-sin2^2)
  refracted <- acos(cos2)  
  refracted
}
#compute angles of refraction  in tissue 2 or reflection
#given incidence angles coming from tissue 1
new_angles <- function(theta,n1,n2){
  angles <- numeric(length(theta))
  reflprobs <- prob_refl(theta,n1,n2)
  idx_refl <- reflprobs==1.0
  sin2 <- n1[!idx_refl]/n2[!idx_refl] * sin(theta[!idx_refl])
  cos2 <- sqrt(1-sin2^2)
  angles[idx_refl] <- theta[idx_refl] 
  angles[!idx_refl] <- acos(cos2)
  list(angle=angles,refl=idx_refl)
}
#given nx3 array of  x,y,z positions
#return nx3 array of voxel (i,j,k) indices
#use ceiling since R is 1-based indexing
get_voxel <- function(P){
#  P[,3] <- P[,3]*.2  #pixels are 1x1x1 instead of 1x1x5
  P <- ceiling(P)  
  P
}
phantomize <- function(x){
  return(as.integer(phantom[x[1],x[2],x[3]]))
}
#this will be used when we consider pixels with dyes
is_stained <- function(x){
  return(as.integer(phantom[x[1],x[2],x[3]])>99)
}
#given nx3 array of voxel indices
#return n array of tissue type
get_tissuetype <- function(V){
  n <- nrow(V)
  tissue <- rep(0,n)
  tissue <- apply(V,1,phantomize)
  tissue
}
#given nx3 array of voxel indices
#return n array of boolean indicating
#which voxels are stained
get_stained <- function(V){
  n <- nrow(V)
  stained <- rep(0,n)
  stained <- apply(V,1,is_stained)
  stained
}
# candir is nx3 array of differences in voxel coordinates
#check if only 1 coord differs
check_canonical <- function(candir){
  #if any rows have more than 1 differing coordinate
  #pick one of these at random and zeroize the rest
  more <- rowSums(candir[1:nrow(candir),]!=0) > 1
  for (i in 1:nrow(candir)) if (more[i])    candir[i,-sample(which(candir[i,]!=0),1)] <- 0
  candir
}
# Simulates scattering and absorption in a phantom head
# assuming nphotons are emitted at the skull with g=.9  
sim_forward <- function(nphotons, g=0.9, max_steps=10){
  steps <- max_steps
  state = list(P = cbind(x=rep(0,nphotons), y=rep(0,nphotons), z=rep(0,nphotons)),
               D = rusphere(nphotons))
  trax  <- list(state$P)
  alive <- matrix(0,nphotons,max_steps)
  dead <- numeric(0,nphotons)
  alive[,1] <- 1
  for(i in 1:steps){
    if(length(dead)==nphotons)break
    nxt <- step(state$P, state$D)
    
    
    temp <- matrify(nxt$X, nxt$X[,"z"] > 0)
    top_exits <- rbind(top_exits, temp[,c("x", "y")])
    state <- nxt[c("P", "D")]
    trax <- c(trax,state$P)
  }
  state$steps <- i
  state
}
tissue_char <- as.data.frame(matrix(c(
  #id & \mu_a & \mu_s   & g &     n &     W   # type\\
  0 ,  0      , 0    ,  .00001 , 1.0    , 0,     #Background
  1 ,  .0076  , 0.01 ,   0.9   , 1.33   , 1.0,   #CSF
  2 ,  0.0335 , 10   , .9      , 1.3688 , .8,    #Grey
  3 ,  0.0207 , 33   , .88     , 1.3852 ,  .7,   #White
  4 ,  .087  ,  0   , 0       , 1.48   ,  0,    #Fat
  5 ,  1.12   ,  53  , .95       , 1.0    ,  0,    #Muscle
  6 , .35     , 35   , .8      , 1.0    , 0,     #Muscle/Skin
  7 ,  .015   , 8.6  , .94      , 1.0    , 0,     #Skull
  8 ,  .5     , 141.3, .99     , 1.0    , 0,     #Vessels
  9 , .087      , 0    , 0       , 1.46    , 0,       #Around Fat (collagen)
  10,  .085   , 6.5  , 0.765   , 1.0    , 0,       #Dura Mater
  11,  .015   , 9.6  , .9      , 1.0    , 0        #Bone Marrow
), 12, 6, byrow=TRUE))
names(tissue_char) <- c("id",  "mu_a","mu_s","g","n","W") 
#z bottom to top
#y back to front
#x left to right
#given 2 nx3 arrays of source and dest points
find_intersects <- function(P1,P2){
  n <- nrow(P1)
  stor <- lapply(1:n,function(x){findVoxelCrossings(P1[x,],P2[x,])})
  stor
}
#Given n-long list of crossing points, 
#n source and n provisional destination points,
#and origin tissue type
#see if crossings go into different tissues
process <- function(xing,src,dest,tissue1){
  n <- length(src)
  temp <- matrix(0,n,3)
  midpt <- matrix(0,1,3)
  idim <- numeric(n)
  ndir <- numeric(n)
  for (i in 1:n)
    if (!is.null(xing[[i]])){
      #compute midpoints between crossings to find voxels correctly
      midpt[1,1:3] <- (src[i,1:3] + xing[[i]][1,1:3])/2
      if (nrow(xing[[i]])>1){
        for (j in 2:nrow(xing[[i]])) rbind(midpt, (xing[[i]][j-1,1:3]+xing[[i]][j,1:3])/2)
        rbind(midpt, (xing[[i]][nrow(xing[[i]]),1:3]+dest[i,1:3])/2)
      } else rbind(midpt, (xing[[i]][1,1:3]+dest[i,1:3])/2)
      #Pi crossed into different voxels, so see what tissues they are
      tissue2<- get_tissuetype(get_voxel(midpt))
      print(paste("photon is ",i," dest tissue is ",tissue2))
      if (sum(tissue1!=tissue2) >0){
        #photon i hits at least 1 different tissue
        #find first voxel crossing with diff tissue
        idx <- which(tissue1!=tissue2)[1]
        print(paste("photon is ",i," idx is ",idx))
        temp[i,1:3] <- xing[[i]][idx,1:3]
        idim[i] <- xing[[i]][idx,4]
        if (temp[i,idim[i]] > src[i,idim[i]]) {
          #aim to greater photon in correct dim
          temp[i,idim[i]] <- temp[i,idim[i]]+1
          ndir[i] <- 1
        }
        else ndir[i] <- -1 #aim back
      } else {#all voxels of same tissue type
        temp[i,1:3] <- dest[i,1:3]
        idim[i] <- 0
        ndir[i]<- 0
      }
    }    else {#no crsossing into diff voxels
      temp[i,1:3] <- dest[i,1:3]
      idim[i] <- 0
      ndir[i]<- 0
    }
  list(P=temp,idim=idim,ndir=ndir)
}
###
# Given:
#   P, an nx3 matrix of internal positions of photons
#   D, an nx3 matrix of unit vectors indicating directions of motion
#   invCDF, the inverse CDF of the scattering angle cosine
# simulate scattering, absorption and exit (contact with boundary) events, returning
# new P and D for remaining photons, and a matrix of exit positions X.
step <- function(P, D){
  n <- nrow(P)
  V <- get_voxel(P)
  # get tissue types of current positions
  tissue1 <- get_tissuetype(V)
  mu_s <- numeric(n)
  mu_a <- numeric(n)
  g <- numeric(n)
  n1 <- numeric(n)

  #get scattering and absorption coeff and index of refraction for source tissues
  #add 1 to tissue1 because tissues are labeled 0..n and lists are indexed 1..n+1
  mu_s <- tissue_char$mu_s[1+tissue1]
  mu_a <- tissue_char$mu_a[1+tissue1]
  n1 <- tissue_char$n[1+tissue1]
  g <- tissue_char$g[1+tissue1]
  invCDF <- lapply(g,icdfHG)

  # step provisionally, ignoring voxel boundaries
  provisional_step <- move_provisionally(P, D, mu_s, mu_a)

  #see what boundaries the provisional steps of photons have crossed
  xing <- find_intersects(P,provisional_step$P)
  
  #change destination points if there were crossings to diff tissues
  newP <- process(xing,P,provisional_step$P,tissue1)
  VP <- get_voxel(newP$P)
  tissue2 <- get_tissuetype(VP)
  n2 <- numeric(length(newP))
  n2 <- tissue_char$n[1+tissue2]
  #check if destination tissue is different from source tissue
  diff_tiss <- tissue1 != tissue2
  #first get canonical direction of movement (left/right, above/below, front/behind)
  candir <- VP[diff_tiss,]-V[diff_tiss,]
  #TODO check candir  to see if there's more than one canonical direction
  candir <- check_canonical(candir)
  src_angles <- acos(rowSums(candir*D[diff_tiss,]))
  #compute new angles of reflection and refraction
  new_ang <- new_angles(src_angles,n1[diff_tiss],n2[diff_tiss])
  #compute directions for all photons(reflected and refracted)
  #which hit different tissue types
  D[diff_tiss] <- compute_direction(new_ang,candir,newP$idim[diff_tiss],provisional_step$D[diff_tiss])
  #check to see if any photons have hit background voxels
  exits <- mark_exits(newp$P,new_ang$refl)
#   # update positions
#   P <- corrected_step[["P"]]
#   # extract exit positions
#   exits <- corrected_step[["exits"]]
   X <- matrify(P, exits)
#   # extract absorbed positions with are not exits
#   absorbed <- provisional_step[["absorptions"]] & !exits
   P <- newP$P
   absorbed <- provisional_step$absorptions & !diff_tiss & !exits
   A <- matrify(P, absorbed)
   # remove absorbed photons from provisional_step components
   P <- matrify(P, !absorbed)
   D <- matrify(D, !absorbed)
  
  # scatter the remaining photons
  D[!absorbed] <- scatter(D[!absorbed], invCDF[!absorbed])

  # return P, D, X, and A
  #list(P=P, D=D, X=X, A=A)

  list(P=P, D=D, A=A)
}
#given list of new angles (and associated boolean indicator of reflection),
#canonical directions and dimensions and source directions compute new directions
#new angles are either reflection or refractions indicated by new_angles$refl
compute_direction <- function(new_angles,candir,wdim,D){
  #form vector of +/- 1's from candir and dimension selector wdim
  vdir <- sapply(1:nrow(candir,function(n){candir[n,wdim[n]]}))
  #multiply angle by +1 or -1 to indicate direction and compute cosines
  newdir <- cos(new_angles$angle*vdir)
  beta <- sqrt(1-newdir^2)
  #create Boolean of refracted angles
  refraction <- !new_angles$refl  
  #0 out canonical directions in D for refracted angles
  for(i in 1:nrow(D)) if (refraction[i]) D[i,wdim[i]] <- 0
  #compute normalizing factor sqrt of sums of sqrs of noncanonical directions
  denom <- sqrt(rowSums(D[refraction,]^2))
  D[refraction,] <- D[refraction,]*beta/denom
  #insert canonical direction
  for (i in 1:nrow(D)) if (refraction[i]) D[i,wdim[i]] <- newdir
  #for reflected angles negate sign of canonical direction
  for (i in 1:nrow(D)) if (!refraction[i]) D[i,wdim[i]] <- -D[i,wdim[i]]
  D
}

# R subsetting casts a 1xm matrix to an m-long vector because, after
# all, consistency is the hobgoblin of small minds. This small-minded
# function subsets consistently, using the logical vector idx to
# return a sub-MATRIX of select rows of M.
matrify <- function(M, idx, ncol=3){
  ans <- matrix(M[idx,], ncol=ncol)
  colnames(ans) <- colnames(M)
  try(rownames(ans) <- rownames(M)[idx], silent=FALSE)
  ans
}
# Given a scattering coefficient (not a reduced scattering coefficient,) mu_s,
# and an absorption coefficient, mu_a, both in units of events per mm, and
# given nx3 arrays, P and D, representing positions and directions of travel
# respectively, compute new positions based on randomly sampled
# distances to new events, assuming these occur within the medium of interest.
# Return the new positions along with a logical vector indicating which events
# were absorptions.
# NOTE: The rows of D must be unit vectors.
move_provisionally <- function(P, D, mu_s, mu_a){
  n <- nrow(P)
  scattering_distances <- rexp(n, mu_s)
  absorption_distances <- rexp(n, mu_a)
  absorptions <- scattering_distances > absorption_distances
  P <- P + pmin(scattering_distances, absorption_distances)*D
  list(P=P, absorptions=absorptions)
}
# Given nx3 matrices of positions, P, and directions of motion, D,
# toward those positions, create a logical vectors marking rows
# for which P is in a background voxel
# Return  exit indicators.
mark_exits <- function(P,refl){
  V <- get_voxel(P)
  tissue <- get_tissuetype(V)
  exits <- (tissue==0)&(!refl)
  exits
}
