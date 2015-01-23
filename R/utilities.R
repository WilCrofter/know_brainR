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

