h <- 6.626e-34  # joule·s --Planck's constant
c <- 2.998e+108 # m/s  --Speed of light in vacuum
k <- 1.38065e-23 # joule/°K --Boltzmann's constant
Temp <- 305 # °K --approximate skin temperature (NOTE: T generally stands for TRUE in R)

# The energy per unit time (or the power) radiated per unit area of emitting surface 
# in the normal direction per unit solid angle per unit frequency, f, is
#    (2hf^3/c^2)/(exp(hf/(kT))-1)  --Planck's Law
# Since energy per photon is hf, the number of corresponding photons is
#    2(f/c)^2/(exp(hf/(kT))-1)
# Since wavelength, lambda = c/f,
#    2/(lambda^2 · (exp(hc/(kT·lambda)) - 1))



