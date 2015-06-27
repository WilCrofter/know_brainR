if(!require("deSolve"))stop("Hodgkin-Huxley code requires package deSolve. Please install it.")

#' In terms applicable to the deSolve package, the Hodgkin-Huxley model has parameters consisting
#' of membrane capacitance, leak conductance (a constant,) and reverse potentials for excitation, 
#' inhibition, and leakage. It has a state consisting of membrane potential alone.
#' Excitation and inhibition conductances are driving functions and, as such, must be incorporated 
#' in a function which returns the derivative of the state variable, membrane potential, as 
#' a function of time, the state, and parameters.
#' 
#' A natural form for a driving function is an 
#' R function with a single argument, time.




#' The default conductance and reverse potential values used below are taken from Destexhe et. al.,
#' FLUCTUATING SYNAPTIC CONDUCTANCES RECREATE INVIVO-LIKE ACTIVITY IN NEOCORTICAL NEURONS
#' http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3320220/.
#' The reference develops a single compartment model for background synaptic activity of a
#' single neuron and shows that it accurately describes the collective behavior of a much more
#' detailed model involving individual synapses. These values thus apply to a single neuron, 
#' and not necessarily as used here in a single compartment model of multiple neurons mixed 
#' with remote and local afferents.

# Conductance units are in micro Siemens per square cm
# (Siemens are amps of current per volt of potential.)
#  Excitation conductance
gE_mean <- 0.012
gE_sd <- 0.003
#  Inhibition conductance
gI_mean <- 0.057
gI_sd <- 0.0066
#  Leak conductance (a constant)
gL <- 0.045

# Membrane capacitance, 1 micro Farad per square sm
C <- 1

# Normal distribution truncated to [0, Inf) by rejection of negative samples.
# This will approximate a formal truncated normal with given mean and variance
# provided the mean is several standard deviations positive.
rtnorm <- function(n, mean=3, sd=1){
  v <- numeric()
  k <- 0
  while(k < n){
    u <- rnorm(n-k, mean, sd)
    u <- u[u >= 0]
    v <- c(v, u)
    k <- length(v)
  }
  return(v)
}

# from SEJ RC


# C m = 1 μF/cm 2 is the specific membrane capacitance, g L = 0.045 mS/cm 2 is the leak
# conductance density, and E L = −80 mV is the leak reversal potential. I Na is the voltage-
#   dependent Na + current and I Kd is the ‘delayed-rectifier’ K + current responsible for action
# potentials. I M is a non-inactivating K + current responsible for spike frequency adaptation.
# These currents and their parameters were the same as in the biophysical model (see
#                                                                                Destexhe and Paré, 1999). a is the total membrane area, which was 34 636 μm 2 for the layer
# VI cell described in Fig. 1.

# E e = 0 mV and E i = −75 mV are their respective reversal potentials and were
# identical to that of the detailed biophysical model.