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
 # print(paste("prob ",prob))
}
#From ../jacques.Rmd
dry <- 1.514
wet <- 1.33
water_gray <- .8
water_white <- .7
water_csf <- 1.0
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
  P[,3] <- P[,3]*.2
  P <- ceiling(P)  
  P
}
phantomize <- function(x){
  return(as.integer(phantom[x[1],x[2],x[3]]))
}
is_stained <- function(x){
  return(as.integer(phantom[x[1],x[2],x[3]])>99)
}
#given nx3 array of voxel indices
#return n array of tissue type
get_tissuetype <- function(V){
  n <- nrow(V)
  tissue <- rep(0,n)
#  for (i in 1:n) tissue[i] <- as.integer(phantom[V[i,1],V[i,2],V[i,3]])
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
#check 2 conditions: only 1 coord differs, and diff tissue voxel is adjacent to source voxel 
# candir is nx3 array of differences in voxel coordinates
check_canon <- function(candir){
  #see if any rows have more than 1 differing coordinate
  #pick one of these at random and zeroize the rest
  more <- rowSums(candir[1:nrow(candir),]!=0) > 1
  pick <- sample(which(candir[more]!=0),1)
  candir[more,-pick] <- 0
  #see if any entries in candir have abs val > 1
  more <- abs(candir[1:nrow(candir),])>1
  if (candir[more]<0) candir[more] <- -1
  else candir[more] <- 1
  candir
}
# Simulates scattering and absorption in a uniform material of infinite extent,
# assuming nphotons are emitted at the skull with g=.9  
sim_forward <- function(nphotons, g=0.9, max_steps=10){
  invCDF <- icdfHG(g)
  steps=max_steps
  state = list(P = cbind(x=rep(0,nphotons), y=rep(0,nphotons), z=rep(0,nphotons)),
               D = rusphere(nphotons))
  for(i in 1:steps){
    if(length(state$P)==0)break
    nxt <- step(state$P, state$D, invCDF)
    temp <- matrify(nxt$X, nxt$X[,"z"] > 0)
    top_exits <- rbind(top_exits, temp[,c("x", "y")])
    state <- nxt[c("P", "D")]
  }
  state$steps <- i
  state
}
tissue_char <- as.data.frame(matrix(c(
  # id & type & \mu_a & \mu_s & g & n & W\\
  0 ,  0 , 0 , 0 , 1.0 , 0,              #Background
  1 ,  .0076 ,0,0, 1.33 , 1.0,           #CSF
  2 ,  0.0335 , 10, .9  , 1.3688 , .8,   #Grey
  3 ,  0.0207 , 33 ,.88 , 1.3852, .7,    #White
  4 ,  .0005 ,0,0, 1.48,0,               #Fat
  5 ,  1.12,53 ,0,1.0,0,                 #Muscle
  6 , .35 , 35, .8 ,1.0,0,               #Muscle/Skin
  7 ,  .015 , 8.6, .9,1.0,0,             #Skull
  8 ,  .5, 141.3, .99,1.0,0,             #Vessels
  9 ,  0, 0, 0, 1.0, 0,                  #Around Fat
  10,  0, 0, 0, 1.0, 0,                  #Dura Mater
  11,  0, 0, 0, 1.0, 0                   #Bone Marrow
), 12, 6, byrow=TRUE))
names(tissue_char) <- c("id",  "mu_a","mu_s","g","n","W") 
#z bottom to top
#y back to front
#x left to right


# Given:
#   P, an nx3 matrix of internal positions of photons
#   D, an nx3 matrix of unit vectors indicating directions of motion
#   invCDF, the inverse CDF of the scattering angle cosine
# simulate scattering, absorption and exit (contact with boundary) events, returning
# new P and D for remaining photons, and a matrix of exit positions X.
step <- function(P, D, invCDF){
  n <- nrow(P)
  # get tissue types of current positions
  V <- get_voxel(P)
  tissue1 <- get_tissuetype(V)
  mu_s <- numeric(n)
  mu_a <- numeric(n)
  n1 <- numeric(n)
  n2 <- numeric(n)
  #add 1 to tissue1 because tissues are labeled 0..n and lists are indexed 1..n+1
  mu_s <- tissue_char$mu_s[1+tissue1]
  mu_a <- tissue_char$mu_a[1+tissue1]
  n1 <- tissue_char$n[1+tissue1]
  # step provisionally, ignoring boundaries
  provisional_step <- move_provisionally(P, D, mu_s, mu_a)
  VP <- get_voxel(provisional_step$P)
  tissue2 <- get_tissuetype(VP)
  n2 <- tissue_char$n[1+tissue2]
  #check for reflection and refraction if destination tissue is different from source tissue
  diff_tiss <- tissue1 != tissue2
  #move P of provisionally stepped photons with different tissues to appropriate boundary
  provisional_step$P[diff_tiss] <-
  #first get canonical direction of movement (left/right, above/below, front/behind)
  candir <- V[diff_tiss,]-VP[diff_tiss,]
  #check candir for 2 conditions
  candir <- check_canonical(candir)

  src_angles <- acos(candir*provisional_step$D[diff_tiss])
  dest_angles <- new_angles(src_angles,n1[diff_tiss],n2[diff_tiss])


  # correct for exits
  corrected_step <- mark_exits(provisional_step[["P"]], D, thickness)
  # update positions
  P <- corrected_step[["P"]]
  # extract exit positions
  exits <- corrected_step[["exits"]]
  X <- matrify(P, exits)
  # extract absorbed positions with are not exits
  absorbed <- provisional_step[["absorptions"]] & !exits
  A <- matrify(P, absorbed)
  # remove exiting and absorbed photons from P and D
  dead <- exits | absorbed
  P <- matrify(P, !dead)
  D <- matrify(D, !dead)
  # scatter the remaining photons
  D <- scatter(D, invCDF)
  # return P, D, X, and A
  list(P=P, D=D, X=X, A=A)
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
# for which the z coordinate exceeds thickness/2 in absolute
# value indicating exit from the tissue. Adjust relevant positions to 
# their exit points. Return top and bottom exit indicators and
# adjusted P.
mark_exits <- function(P, D, thickness){
  tops <- P[,"z"] > thickness/2
  bottoms <- P[ ,"z"] < -thickness/2
  P[tops,] <- P[tops, ] - ((P[tops, "z"]-thickness/2)/D[tops, "z"])*D[tops,]
  P[bottoms,] <- P[bottoms, ] - ((P[bottoms, "z"]+thickness/2)/D[bottoms, "z"])*D[bottoms,]
  list(P=P, exits=tops | bottoms)
}
