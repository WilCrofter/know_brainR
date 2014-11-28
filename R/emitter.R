#' emit n photons from a given spectrum
#' 
#' @param n number of photons, a positive integer
#' @param spectrum a data frame or matrix whose first column is wavelength 
#' and whose second column is the probability of emission at the corresponding wavelength.
#' @return a vector of n emission wavelengths
emit <- function(n, spectrum){
  sample(as.vector(spectrum[,1]), n, replace=TRUE, prob=as.vector(spectrum[,2]))
}