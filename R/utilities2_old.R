#myseed <- 0x1257AB #run300
#myseed <- 0x7141958 #lp301
#myseed <- 0x8728672 #run500
#myseed <- 0x2487ad #run501
set.seed(myseed)

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
#TODO put in concentration of dyes
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
  more <- rowSums(abs(candir) > 0)
  for (i in 1:nrow(candir)) if (more[i]>1)  candir[i,-sample(which(candir[i,]!=0),1)] <- 0
  candir
}
# run from directory know_braimR
read_phantom <- function(){
fname <- "data/subject04_crisp_v.rawb"
# Read in raw bytes as a vector
phantom <- readBin(fname, what="raw", n=362*434*362, size=1, signed=FALSE, endian="big")
# Convert to a 3D array by setting the dim attribute
dim(phantom) <- c(362, 434, 362)
phantom
}


# Simulates scattering and absorption in a phantom head
# assuming nphotons are emitted at the skull with g=.94 
sim_forward <- function(nphotons, max_steps=1000){
  steps <- max_steps
  print(paste("seed is",myseed))

  #state holds position, direction,and condition, i.e., alive, absorbed,exit of each photon
  state = list(
    P = cbind(x=rep(109.5,nphotons), y=rep(216.5,nphotons), z=rep(327.5,nphotons)),
    D = cbind(x=rep(-4/45.177,nphotons), y=rep(0,nphotons), z=rep(-45/45.177,nphotons)), 
    flg=rep("Alive",nphotons))
#    P = cbind(x=rep(181,nphotons), y=rep(217,nphotons), z=rep(340,nphotons)),
#     P = cbind(x=rep(180.5,nphotons), y=rep(216.5,nphotons), z=rep(339.5,nphotons)),
#     D = cbind(x=rep(0,nphotons), y=rep(0,nphotons), z=rep(-1,nphotons)), flg=rep("Alive",nphotons))
  trax  <- list(cbind(state$P, state$flg))
  out <- list(cbind(x=numeric(),y=numeric(),z=numeric()))
  absorbed <- list(cbind(x=numeric(),y=numeric(),z=numeric()))

  record <- matrix("-",nphotons,max_steps+1)
  path <- cbind(state$P, state$flg)

  record[,1] <- "Alive"
  for(i in 1:steps){
    if(sum(record[,i]=="Alive")==0)break
#    print(paste("step = ",i,"num alive is ",sum(record[,i]=="Alive")))
    nxt <- steps(state$P, state$D, state$flg)
    
    out[[i+1]] <- nxt$X
    absorbed[[i+1]] <- nxt$A
    record[record[,i]=="Alive",i+1] <- nxt$flg
  
    path[record[,i]=="Alive",] <- c(nxt$P[,1:3],nxt$flg[])
    path[record[,i]!="Alive",] <- c("-","-","-","-")

    nxt$P <- matrify(nxt$P,nxt$flg=="Alive")
    nxt$D <- matrify(nxt$D,nxt$flg=="Alive")
    
 
    
    nxt$flg <-nxt$flg[nxt$flg=="Alive"]    
    state <- nxt[c("P", "D", "flg")]
    trax[[i+1]] <- path
  }
#  state
#  print(record)

list(trax=trax,record=record,tchar=tissue_char,seed=myseed)


}
num_types <-12
num_char <- 6
get_tissue_chars <- function(){
  means <- matrix(c(
  #id & \mu_a & \mu_s   & g &     n &     W   # type\\
  0 ,  0.00001, 0.00001, .00001 , 1.0    , 0,     #Background
  1 ,  .0076  , 0.01 ,   0.9   , 1.33   , 1.0,   #CSF
  2 ,  0.0335 , 10   , .9      , 1.3688 , .8,    #Grey
  3 ,  0.0207 , 33   , .88     , 1.3852 ,  .7,   #White
  4 ,  .087  ,  11.5   , .9       , 1.48   ,  0,    #Fat
  5 ,  1.12   ,  53  , .95       , 1.41    ,  0,    #Muscle
  6 , .35     , 35   , .93      , 1.45    , 0,     #Muscle/Skin g=.8?
  7 ,  .015   , 8.6  , .94      , 1.55    , 0,     #Skull
  8 ,  .22     , 58.5, .99     , 1.4    , 0,     #Vessels
  9 , .087      ,11.5    , .9       , 1.46    , 0,       #Around Fat (collagen)
  10,  .085   , 6.5  , 0.765   , 1.4    , 0,       #Dura Mater
  11,  .015   , 9.6  , .9      , 1.4    , 0        #Bone Marrow
  ), num_types, num_char, byrow=TRUE)

  std_dev <- matrix(0.000001,nrow=num_types,ncol=num_char)
  std_dev[c(9,11,12),5] <- .01
  tissue_char <- as.data.frame(gen_tissue_chars(means,std_dev) )
  names(tissue_char) <- c("id",  "mu_a","mu_s","g","n","W") 
  tissue_char
}#end of get_tissue_chars
gen_tissue_chars <- function(means,std_dev){
  temp <- matrix(0,num_types,num_char)
  for (i in 1:num_types){
    for (j in 2:num_char) temp[i,j] <- rnorm(1,means[i,j],std_dev[i,j])
  }
  temp
}

phantom <- read_phantom()
tissue_char <- get_tissue_chars()

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
  n <- nrow(src)
  temp <- matrix(0,n,3)
  idim <- numeric(n)
  ndir <- numeric(n)
  for (i in 1:n)
    if (!is.null(xing[[i]])){
      #compute midpoints between crossings to find voxels correctly
      midpt <- matrix(0,1,3)
      midpt[1,1:3] <- (src[i,1:3] + xing[[i]][1,1:3])/2
      if (nrow(xing[[i]])>1){
        for (j in 2:nrow(xing[[i]])) rbind(midpt, (xing[[i]][j-1,1:3]+xing[[i]][j,1:3])/2)
        midpt <- rbind(midpt, (xing[[i]][nrow(xing[[i]]),1:3]+dest[i,1:3])/2)
      } else  midpt <- rbind(midpt, (xing[[i]][1,1:3]+dest[i,1:3])/2)
      #Pi crossed into different voxels, so see what tissues they are
      tissue2<- get_tissuetype(get_voxel(midpt))
      if (sum(tissue1[[i]]!=tissue2) >0){
        #photon i hits at least 1 different tissue
        #find first voxel crossing with diff tissue
        idx <- which(tissue1[[i]]!=tissue2)[1]
        temp[i,1:3] <- xing[[i]][idx-1,1:3]
        idim[i] <- xing[[i]][idx-1,4]
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
# flg, an n long character vector indicating "Alive", "Exited" or "Absorbed"
# simulate scattering, absorption and exit (contact with boundary) events, returning
# new P and D for remaining photons, and  matrices of exit positions  and absorbed photons.
steps <- function(P, D, flg){
  n <- nrow(P)
  V <- get_voxel(P)
  # get tissue types of current positions
  tissue1 <- get_tissuetype(V)
    
  mu_s <- numeric(n)
  mu_a <- numeric(n)
  g <- numeric(n)
  n1 <- numeric(n)
  src_angles <- numeric()
  new_ang <- list()

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
  n2 <- numeric(length(newP$idim))
  n2 <- tissue_char$n[1+tissue2]
  #check if destination tissue is different from source tissue
  diff_tiss <- tissue1 != tissue2
  #first get canonical direction of movement (left/right, above/below, front/behind)
  candir <- VP[diff_tiss,]-V[diff_tiss,]
  if (sum(diff_tiss)==1){
     dim(candir) <- c(1,3)
   }
  if ((nrow(candir)>0)){
    #check candir  to see if there's more than one canonical direction in any row
    candir <- check_canonical(candir)
    src_angles[diff_tiss] <- acos(rowSums(candir*D[diff_tiss,]))
    #compute new angles of reflection and refraction
    temp <- new_angles(src_angles[diff_tiss],n1[diff_tiss],n2[diff_tiss])

    #compute directions for all photons(reflected and refracted)
    #which hit different tissue types
    D[diff_tiss,] <- compute_direction(temp,candir,newP$idim[diff_tiss],D[diff_tiss,1:3])
    new_ang <- list(angle=numeric(n),refl=logical(n))
    new_ang$angle[diff_tiss]<-temp$angle
    new_ang$refl[diff_tiss] <- temp$refl
 
  }
  else{
    new_ang <- list(angle=numeric(n),refl=logical(n))
  }
  new_ang$refl[!diff_tiss] <- FALSE

  #if reflected (hence different  tissues), put photon back in direction of source voxel
  #note that dimension that's being changed is an integer, i.e., voxel boundary
  newP$P[new_ang$refl,newP$idim[diff_tiss]] <- newP$P[new_ang$refl,newP$idim[diff_tiss]]-newP$ndir[diff_tiss]
  P <- newP$P
  
  #check to see if any photons have hit background voxels and aren't reflected back
  exits <- mark_exits(P,new_ang$refl)
  flg[exits] <- "Exited"
   # extract exit positions
   X <- matrify(P, exits)
   # extract absorbed positions with are not exits and not different tissues
   absorbed <- provisional_step$absorptions & !diff_tiss & !exits
   flg[absorbed] <- "Absorbed"
   A <- matrify(P, absorbed)
   #mark photons that haven't exited or been absorbed
   alive <- !(exits | absorbed)
   flg[alive] <- "Alive"

   # scatter the remaining photons
   for (i in 1:n) if (flg[i]=="Alive") D[i,1:3] <-scatter1(D[i,1:3],invCDF[[i]](runif(1))) 

   # return P, D, X, A, and indicator flag
   list(P=P, D=D, X=X, A=A, flg=flg)

}
#given list of new angles (and associated boolean indicator of reflection),
#canonical directions and dimensions and source directions compute new directions
#new angles are either reflection or refractions indicated by new_angles$refl
compute_direction <- function(new_angles,candir,wdim,D){
  #check to make sure D is matrix
  if (length(D)==3){
    dim(D) <- c(1,3)
  }
  #form vector of +/- 1's from candir and dimension selector wdim
  vdir <- sapply(1:nrow(candir),function(n){candir[n,wdim[n]]})

  #multiply angle by +1 or -1 to indicate direction and compute cosines
  newdir <- cos(new_angles$angle*vdir)
  beta <- sqrt(1-newdir^2)
  #create Boolean of refracted angles
  refraction <- !new_angles$refl  
  #0 out canonical directions in D for refracted angles
  for(i in 1:nrow(D)) if (refraction[i]) D[i,wdim[i]] <- 0
  #compute normalizing factor sqrt of sums of sqrs of noncanonical directions
  denom <- sqrt(rowSums(D^2))[refraction]
  D[refraction,] <- D[refraction,]*beta/denom
  #insert canonical direction
  for (i in 1:nrow(D)) if (refraction[i]) D[i,wdim[i]] <- newdir[i]
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
#Given mxn array  records (m=number of photons, n=number of steps+1)
#and n-long list where each entry in list is mx4 array
#return an mx5 array of last positions before photon died (absorbed or exited skull)
#first 3 columns indicate position, 4th cause of death if any, 5th time of death
last_position <- function(record,lpos){
  nsteps <- length(lpos)
  nphotons <- nrow(record)

  book <- matrix(0,nphotons,5)
  alive <- which(record[,nsteps]=="Alive")
  book[alive,1:3] <- lpos[[nsteps]][alive,1:3]
  book[alive,4] <- "Alive"
  book[alive,5] <- nsteps 
  for (i in 1:nsteps-1) {
     idx <- which(record[,nsteps-i]!="Alive" & record[,nsteps-i]!="-")
     book[idx,1:4] <- lpos[[nsteps-i]][idx,1:4]
     book[idx,5] <- nsteps-i
  }
  
  book
}
#Given n-long list where each entry in list is mx4 array of photon positions and cause
# of death calculate length of each photon's path
path_length <- function(lpos){
  nsteps  <- length(lpos)
  nphotons <- nrow(lpos[[1]])
  plens <- matrix(0,nphotons,nsteps)
  temp <- rep(0,nphotons)
  pos1 <- matrix(as.numeric(lpos[[1]][,1:3]),nphotons,3)
  pos2 <- matrix(0,nphotons,3)
  for (i in 2:nsteps){
    idx <- lpos[[i]][,1]!="-"  #photons still alive at step i
    len <- sum(idx)
    pos2[idx,1:3] <- matrix(as.numeric(lpos[[i]][idx,1:3]),len,3)
    pos2[!idx,1:3] <- 0
    temp <- rep(0,nphotons)
    temp[idx] <- (pos2[idx,1]-pos1[idx,1])^2 +
      (pos2[idx,2]-pos1[idx,2])^2 +
      (pos2[idx,3]-pos1[idx,3])^2
    plens[idx,i] <- plens[idx,i-1] + sqrt(temp[idx])
    plens[!idx,i] <- plens[!idx,i-1]
    pos1 <- pos2

  }
plens
}
#Given nx5 array (output of last_position) where n=number of photons, 
#cols 1-3 are last position of photon, col 4= cause of death, col 5=time of death
#and P0=starting position of photons, return xy distance from starting position, z distance
#or depth, and time
distance_asfunc <- function(last,P0){
  n <- nrow(last)
  ilast <- matrix(0,n,3)
#  dist <- matrix(0,n,3)
  dist <- cbind(dst=rep(0,n),depth=rep(0,n),step=rep(0,n))
  ilast[,1:3] <- as.numeric(last[,1:3])
  dist[,1] <- sqrt((ilast[,1]-P0[1])^2 + (ilast[,2]-P0[2])^2)
  dist[,2] <- abs(ilast[,3]-P0[3])
  #dist[,3] <- as.integer(last[,5])
  dist
}