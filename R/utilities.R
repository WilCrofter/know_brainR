vsd_spectrum <- function(n, lambda_peak=900, lambda_range=250){
  # Template spectral data
  template <- as.data.frame(matrix(c(
    # x, F, dF/F
    125, 0, 0,
    200, 1.818182e-01, 5.000000e-01,
    240, 5.151515e-01, 6.666667e-01,
    315, 1.000000e+00, 0,
    362, 6.969697e-01, -3.333333e-01,
    440, 3.333333e-01, -2.222222e-01,
    550, 6.060606e-02, -5.555556e-02,
    660, 0, 0
  ), 8, 3, byrow=TRUE))
  names(template) <- c("x", "F", "dlnF") 
  # The peak wavelength corresponds to x=315.
  # The spectrum's range corresponds to 125 < x < 660
  # Therefore the line which converts x to lambda has
  # slope lambda_range/(660-125) and goes through (315, peak_wavelength)
  slope <- lambda_range/(660-125)
  # Define 2 points on the line converting x to wavelength
  xpts <- c(315, 660)
  lpts <- c(lambda_peak, lambda_peak+slope*(660-315))
  # Line defined by the two points
  x2lambda <- lm(lpts ~ xpts)
  # Wavelengths corresponding to x.
  template[,"lambda"] <- x2lambda$coef[1] + x2lambda$coef[2]*template[,"x"]
  # Divide the range of wavelengths into the requested number, n
  rng <- range(template[,"lambda"])
  ans <- data.frame(lambda=seq(rng[1], rng[2], length.out=n))
  # Smoothly interpolate from the template to this range
  ans[,"F_polarized"] <- predict(loess(F ~ lambda, template), ans)
  dlnF <- predict(loess(dlnF ~ lambda, template), ans)
  ans[,"F_depolarized"] <- ans[,"F_polarized"]*(1+dlnF)
  # Zeroize small negatives
  ans[ans < 0] <- 0
  # Normalize intensity curves to describe probabilities of emission.
  ans[,"F_polarized"] <- ans[,"F_polarized"]/sum(ans[,"F_polarized"])
  ans[,"F_depolarized"] <- ans[,"F_depolarized"]/sum(ans[,"F_depolarized"])
  return(ans)
}

# Generates n points uniformly distributed on the unit sphere
# or half sphere (half=TRUE).
# Reference: http://mathworld.wolfram.com/SpherePointPicking.html
rusphere <- function(n, half=FALSE){
  theta <- runif(n, 0, 2*pi)
  cosu <- runif(n, half-1, 1) # 'half' is automatically cast from logical to numerical
  sinu <- sqrt(1-cosu^2)
  cbind(x=sinu*cos(theta), y=sinu*sin(theta), z=cosu)
}

# Henyey-Greenstein approximation to probability density function
# of the cosine of a scattering angle.
pdfHG <- function(cosu, g){
  0.5*(1-g^2)/(1+g*(g-2*cosu))^1.5
}

# Creates cumulative distribution function for the Henyey-Greenstain PDF.
# NOTE: the returned function retains a reference to the lexical scope in
# which it was created, so g and the temporary function will be implicitly
# available. e.g.  if f <- cdfHG(.9), environment(f)$g will return .9.
cdfHG <- function(g){
  temp <- function(cosu)pdfHG(cosu, g)
  function(cosv)sapply(cosv, FUN=function(x)integrate(temp, -1, x)$value)
}

# Creates an inverse CDF function for the Henyey-Greenstein PDF
# associated with anisotropy coefficient g, to enable random
# sampling from the distribution.
icdfHG <- function(g){
  cdf <- cdfHG(g)
  cosu <- seq(-1, 1, by=0.01)
  u <- cdf(cosu)
  approxfun(u, cosu)
}

# Given a unit vector, d, representing direction of travel,
# and the cosine of a scattering angle, cosu, returns a new unit
# vector representing direction after scattering.
scatter1 <- function(d, cosu){
  sinu <- sqrt(1 - cosu^2)
  # orthonormal basis for the orthogonal complement of d
  b1 <- c(d[2], -d[1], 0)/sqrt(sum(d[1:2]^2))
  b2 <- c(0, -d[3], d[2]) + b1*b1[2]*d[3]
  b2 <- b2/sqrt(sum(b2^2))
  # random angle
  psi <- runif(1, 0, 2*pi)
  cosu*d + sinu*sin(psi)*b1 + sinu*cos(psi)*b2
}

# Given an nx3 matrix, D, of unit vectors representing directions
# of photons traveling within a medium with an inverse CDF, invCDF,
# for scattering angles, compute scattered directions.
scatter <- function(D, invCDF){
  # randomly sample from the scattering angle distribution
  # and append cosines of scattering angles to D
  D <- cbind(D, invCDF(runif(nrow(D))))
  # apply scatter1 to each row of D and return the result
  t(apply(D, 1, function(d){scatter1(d[1:3], d[4])}))
}

# Given a scattering coefficient (not a reduced scattering coefficient,) mu_s,
# and an absorption coefficient, mu_a, both in units of events per mm, and
# given nx3 arrays, P and D, representing positions and directions of travel
# respectively, compute new positions based on randomly sampled
# distances to new events, assuming these occur within the medium of interest.
# Return the new positions along with a logical vector indicating which events
# were absorptions.
# NOTE: The rows of D must be unit vectors.
move_tentatively <- function(P, D, mu_s, mu_a){
  n <- nrow(P)
  scattering_distances <- rexp(n, mu_s)
  absorption_distances <- rexp(n, mu_a)
  absorptions <- scattering_distances > absorption_distances
  P <- P + pmin(scattering_distances, absorption_distances)*D
  list(positions=P, absorptions=absorptions)
}

