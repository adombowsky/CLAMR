cluster.gibbs.topiclist <- function(R, B, seed,
                                x, miss, 
                                var.priors, 
                                mu.old.center,mu.old.scale,
                                K, L, 
                                alpha_lambda, rho,
                                stops) {
  # libraries
  set.seed(seed)
  library(Rcpp)
  library(RcppArmadillo)
  print("Sourcing C++ Functions, Stand By")
  sourceCpp("rcppfuncts/rdirichlet_arma.cpp")
  sourceCpp("rcppfuncts/cluster_labels.cpp")
  # preliminaries
  n <- nrow(x)
  d <- ncol(x)
  # initialization
  #c <- sample(L, n, replace=TRUE, prob=rep(1/L,L)) # completely random
  c <- kmeans(x, centers=L)$cluster # with k-means
  s <- matrix(0, nrow=L, ncol = d )
  for (j in 1:d) {
    s[,j] <- sample(K[j], size = L, replace = TRUE)
  }
  ll.out <- c()
  c.out <- matrix(0, nrow=R, ncol=n)
  s.out <- array(0, dim = c(R, L, d))
  c.out[1,] <- c
  s.out[1,,] <- s
  xhat <- x
  mu <- sigma_sq <- matrix(0, nrow = L, ncol = d)
  for (l in 1:L) {
    for (j in 1:d) {
      sigma_sq[l,j] <- 1/rgamma(n=1, shape=var.priors[j,1], rate=var.priors[j,2])
      mu[l,j] <- mu.old.center[[j]][s[l,j]] + sqrt(mu.old.scale[[j]][s[l,j]]) * rnorm(1)
    }
  }
  ldens.out <- matrix(0, nrow=R, ncol=n)
  ldens <- c()
  # configuring stops
  stops <- (1:(R/stops)) * stops
  # sampling
  print("Sampling!")
  for (r in 2:R) {
    # sample c
    c <- update_c_marginal(y=x, c=matrix(c), 
                           mu = mu, 
                           sigma_sq = sigma_sq, 
                           alpha_lambda = alpha_lambda)
    for (j in 1:d) {
      s[,j] <- update_c_marginal(y=matrix(mu[,j]), c = s[,j],
                                 mu = matrix(mu.old.center[[j]]),
                                 sigma_sq = matrix(mu.old.scale[[j]]),
                                 alpha_lambda = rho[j])
    }
    # sample mean and variance
    for (l in 1:L) {
      n.prime.l <- sum(c==l)
      if (n.prime.l == 0) { # if cluster empty, we generate sample from prior
        for (j in 1:d) {
          sigma_sq[l,j] <- 1/rgamma(n=1, shape=var.priors[j,1], rate=var.priors[j,2])
          mu[l,j] <- mu.old.center[[j]][s[l,j]] + sqrt(mu.old.scale[[j]][s[l,j]]) * rnorm(1)
        }
      } else if (n.prime.l == 1) { # if singleton cluster, have to be careful in R
        x.prime.l <- x[c==l, ]
        for (j in 1:d) {
          # input prior parameters
          mu0 <- mu.old.center[[j]][s[l,j]]
          tau0 <- mu.old.scale[[j]][s[l,j]]
          alpha <- var.priors[j,1]
          beta <- var.priors[j,2]
          x.bar.l.j <- x.prime.l[j]
          # first, update sigma^2
          alpha.star <- alpha + 0.5 * n.prime.l
          beta.star <- beta + 0.5 * (x.prime.l[j] - mu[l,j])^2
          sigma_sq[l,j] <- 1/rgamma(n=1, shape=alpha.star, rate = beta.star)
          # now update mu
          tau0.star <- 1/( (1/tau0) + (n.prime.l/sigma_sq[l,j]))
          mu0.star <- tau0.star * ( (mu0/tau0)  + (n.prime.l * x.bar.l.j/sigma_sq[l,j]))
          mu[l,j] <- mu0.star + sqrt(tau0.star) * rnorm(1)
        }
      } 
      else {
        x.prime.l <- x[c==l, ]
        for (j in 1:d) {
          # input prior parameters
          mu0 <- mu.old.center[[j]][s[l,j]]
          tau0 <- mu.old.scale[[j]][s[l,j]]
          alpha <- var.priors[j,1]
          beta <- var.priors[j,2]
          x.bar.l.j <- mean(x.prime.l[,j])
          # first, update sigma^2
          alpha.star <- alpha + 0.5 * n.prime.l
          beta.star <- beta + 0.5 * sum((x.prime.l[,j] - mu[l,j])^2)
          sigma_sq[l,j] <- 1/rgamma(n=1, shape=alpha.star, rate = beta.star)
          # now update mu
          tau0.star <- 1/( (1/tau0) + (n.prime.l/sigma_sq[l,j]))
          mu0.star <- tau0.star * ( (mu0/tau0)  + (n.prime.l * x.bar.l.j/sigma_sq[l,j]))
          mu[l,j] <- mu0.star + sqrt(tau0.star) * rnorm(1)
        }
      }
    }
    
    ## finally, impute missing data
    for (i in 1:n) {
      if (sum(miss[i,] == 0)) {
        next
      } else{
        for (j in 1:d) {
          if (miss[i,j]) {
            x[i,j] <- mu[c[i], j] + sqrt(sigma_sq[c[i], j]) * rnorm(1)
          }
        }
      }
    }
    # sample log-likelihood (for WAIC calculations)
    ll <- 0
    ldens <- rep(0,n)
    for (i in 1:n) {
      for (j in 1:d) {
        ll <- ll + dnorm(x[i,j],mu[c[i],j], sqrt(sigma_sq[c[i],j]), log=T)
        ldens[i] <- ldens[i] + dnorm(x[i,j],mu[c[i],j], sqrt(sigma_sq[c[i],j]), log=T)
      }
    } 
    
    
    # save output
    c.out[r, ] <- as.vector(c)
    s.out[r, ,] <- s
    ll.out[r] <- ll
    ldens.out[r,] <- ldens
    xhat <- xhat + x
    
    # print stops
    if (r %in% stops){
      print(r)
    }
  }
  xhat <- xhat/R
  return(list(c=c.out[-(1:B),], xhat = xhat, s = s.out[-(1:B), ,], ll = ll.out[-(1:B)], ldens = ldens.out[-(1:B),]))
}