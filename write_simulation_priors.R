##### new prior elicitation #####
write_sim_priors <- function(p, K, l_MRs, u_MRs,xi0,tau0) {
  B <- matrix(0, nrow = p, ncol = 2)
  for (j in 1:p) {
    B[j,] <- c(l_MRs[1,j], u_MRs[K,j])
  }
  mu.old.center <- mu.old.scale <- list()
  for (j in 1:p) {
    mu.old.center[[j]] <- mu.old.scale[[j]] <- rep(0, K)
  }
  for (j in 1:p) {
    mu.old.center[[j]] <- xi0[,j]
    mu.old.scale[[j]] <- tau0[,j]
  }
  
  ##### and we derive the variance priors #####
  var.priors <- matrix(0, nrow = p, ncol = 2)
  H <- 10
  range.vec <- (B[,2] - B[,1])/H
  inv.range.vec <- 1/range.vec
  var.priors[,1] <- (inv.range.vec) * H
  var.priors[,2] <- H
  
  ##### create mu priors #####
  mu.priors <- t(apply(B, 1, function(y) interval_to_gaussian(y[1],y[2]) )) 
  return(list(mu.old.center = mu.old.center, mu.old.scale = mu.old.scale, var.priors = var.priors))
}

