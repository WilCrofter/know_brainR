#This is sample code taken from
#http://blog.abhranil.net/2014/02/08/r-code-for-multivariate-random-walk-metropolis-hastings-sampling/

# Define arguments
ring2D <- function(x)    # x is a vector
{
  exp(-5*abs(x[1]^2+x[2]^2-1))
}

vcov2D <- .01*diag(2)

ring3D <- function(x)    # x is a vector
{
  exp(-5*abs(x[1]^2+x[2]^2+x[3]^2-1))
}

#variance-covariance matrix is, say, .01 down the diagonal again:
  vcov3D <- .01*diag(3)


#Define the sampling function
rwmetro <- function(target,N,x,VCOV,burnin=0)
{
  require(MASS)   #requires package MASS for normal sampling
  samples <- x
  for (i in 2:(burnin+N))
  {  
    prop <- mvrnorm(n = 1, x, VCOV)
    if (runif(1) < min(1, target(prop)/target(x)))
      x <- prop
    samples <- rbind(samples,x)
  }
  samples[(burnin+1):(N+burnin),]
}

# Call the sampling function with the arguments
ringsample2D <- rwmetro(ring2D,40000,c(0,0),vcov2D)

# Use the sample
plot(ringsample2D[,1], ringsample2D[,2], xlim=c(-1.5,1.5),ylim=c(-1.5,1.5), main="Metropolis-Hastings Sample",xlab="x", ylab="y", pch='.')

#Let’s go for a sample of 20000 points after a burn-in of 20000, 
#starting from (0,0,0). The function is then called as:

#ringsample3D <- rwmetro(ring3D,20000,c(0,0,0),vcov3D,20000)

#ringsample3D is now a 3×1000 matrix, and I use the following to get a 
#rotatable 3D plot of it (requires the rgl package):
#  require(rgl)

#plot3d(ringsample3D[,1], ringsample3D[,2], ringsample3D[,3], xlab='x',ylab='y', zlab='z',col="red", size=1)
