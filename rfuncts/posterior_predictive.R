posterior_predictive <- function(c.samps, mu.samps,sigma_sq.samps,stops) {
  R <- nrow(c.samps) # number of MCMC iterations
  n <- ncol(c.samps) # sample size
  p <- dim(mu.samps)[2] # number of features
  x_out <- array(0,dim=c(R,n,p))
  x <- matrix(0,nrow=n,ncol=p)
  for (r in 1:R) {
    for (i in 1:n) {
      for (j in 1:p) {
        x[i,j] <-  mu.samps[,,r][c.samps[r,i], j] + 
          sqrt(sigma_sq.samps[,,r][c.samps[r,i], j]) * rnorm(1)
      }
    }
    x_out[r,,] <- x
    if (r%%stops==0) {
      paste("Iter",r)
    }
  }
  return(list(x=x_out))
}