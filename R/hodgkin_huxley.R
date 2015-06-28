if(!require("deSolve"))stop("Hodgkin-Huxley code requires package deSolve. Please install it.")

#' In terms applicable to the deSolve package, a single-compartment, Hodgkin-Huxley model has 
#' parameters consisting of membrane capacitance, C, leak conductance, gl (a constant,) and reverse 
#' potentials, Ee, Ei, El, for excitation, inhibition, and leakage. It has a state consisting of 
#' membrane potential, V, alone.
#' 
#' Excitation and inhibition conductances, ge and gi, are thought of as driving forces which
#' consist of ion channels opening in response to signals which originate outside of the compartment
#' and independently of it. For purposes of using the deSolve package, ge and gi must be 
#' incorporated in an R function which returns the derivative of the state variable as a function 
#' of time, state, and parameters. 
#' 
#' A reasonable way to handle these details in R is to create solvers which internalize
#' model parameters and accept ge, gi and time as arguments. The following function
#' returns such a solver.

#' Create a solver for a Hodgkin-Huxley model with the given parameters.
#' Parameter default values are taken from Destexhe et. al.,
#' http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3320220/.
#' @param V0 the initial membrane potential (millivolts)
#' @param C membrane capacitance (microfarads per square cm)
#' @param gl leakage conductance (milliSiemens per square cm)
#' @param Ee excitation reverse potential (millivolts)
#' @param Ei inhibition reverse potential (millivolts)
#' @param El leakage reverse potential (millivolts)
#' @return a function of ge, gi, and t where ge and gi are themselves R functions
#' which return conductance in mS/cm^2 as a function of time and t is a vector
#' of times at which the response is desired.
HH_model <- function(V0=-80, C=1, gl=0.045, Ee=0, Ei=-75, El=-80){
  function(t, ge, gi){
    parameters=c(C=C, gl=gl, Ee=Ee, Ei=Ei, El=El)
    dV=function(t, V, parameters){
      list((-1/C)*(ge(t)*(V-Ee) + gi(t)*(V-Ei) + gl*(V-El)))
    }
    ode(V0, t, dV, parameters)
  }
} 



