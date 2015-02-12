#computes probability of reflection
#given angle of incidence (in radians)
#and two indices of refraction, 
#source and dest media respectively
prob_refl <- function(theta1,n1,n2){
  prob <- rep(1,length(theta1))
  cos2 <- rep(0,length(theta1))
  theta2 <- rep(-1,length(theta1))
  sin1 <- sin(theta1)
  cos1 <- sqrt(1-sin1*sin1)
  sin2 <- (n1/n2) * sin1
  idx <- sin2<=1
  cos2[idx] <- sqrt(1-sin2[idx]^2)
  theta2[idx] <- arccos(cos2[idx])
  t1 <-  (abs((n1[idx]*cos1[idx]-n2[idx]*cos2[idx])/(n1[idx]*cos1[idx]+n2[idx]*cos2[idx])))^2
  t2 <-  (abs((n1[idx]*cos2[idx]-n2[idx]*cos1[idx])/(n1[idx]*cos2[idx]+n2[idx]*cos1[idx])))^2  
  prob[idx] <- pmin(1, .5 * (t1+t2))
  prob
 # print(paste("prob ",prob))
}

Rgraytowhite <- function(theta){
  dry <- 1.14
  wet <- 1.33

  n_gray <- dry - (dry-wet)*rep(.8,length(theta))
  n_white <- dry - (dry-wet)*rep(.7,length(theta) )
  prob_refl(theta,n_gray,n_white)
  
}
#represent tissue with 5 parameters: depth, mu_a, mu_s, g, n