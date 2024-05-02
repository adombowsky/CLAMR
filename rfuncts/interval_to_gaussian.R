interval_to_gaussian <- function(a, b) {
  mu = (a+b)/2
  sigma_sq = ((b-mu)/1.96)^2
  return(c(mu, sigma_sq))
}

interval_to_gaussian_general <- function(a,b,p) {
  mu = (a+b)/2
  sigma_sq = ((b-a)/(2*qnorm((1+p)/2)))^2
  return(c(mu, sigma_sq))
}