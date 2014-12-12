# Structural analysis of the thickness of human dura mater with scanning electron microscopy
# http://www.ncbi.nlm.nih.gov/pubmed/8815466
#
# Calvarial thickness and its relation to cranial bone harvest.
# http://www.ncbi.nlm.nih.gov/pubmed/16651971
# 
# Behaviour of near-infrared light in the adult human head:
# implications for clinical near-infrared spectroscopy
# http://www.ncbi.nlm.nih.gov/pubmed/10740545
#
# Cranial bone 6.3 mm
# Dura .27 mm
# Cerebral Cortex .29 mm

surviving_photons <- function(m, n, b){
  d <- sqrt(m^2 + n^2 + b^2)
  0.92^d*b/(4*pi*d^1.5)
}

add_contribution <- function(scalp, m, n, b){
  for(h in 1:nrow(scalp)){
    for(w in 1:ncol(scalp)){
      scalp[h,w] <- scalp[h,w] + surviving_photons(m-h, n-w, b)
    }
  }
  scalp
}

scalp <- matrix(0, 51, 51)

for(y in 1:51){
  for(x in 1:10){
    scalp <- add_contribution(scalp, x, y, 10)
  }
}

scalp <- scalp/max(scalp)

image(1:51, 1:51, scalp)
