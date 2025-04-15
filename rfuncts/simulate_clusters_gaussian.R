simulate_clusters_gaussian <- function(n, seed) {
  source("rfuncts/interval_to_gaussian.R")
  source("rfuncts/vanilla_gmm.R")
  source("write_simulation_priors.R")
  source("rfuncts/marginal_cluster_list.R")
  #library(sn)
  library(mclust)
  # simulation case
  set.seed(seed)
  p <- 6 # dimension
  L0 <- 3 # true number of clusters
  c0 <- sample(L0, size = n, replace = TRUE) # true clustering
  K <- 3
  ## MRs ##
  l_MRs <- matrix(c(
    -1,0,-10,-0.1, 0, 0,
    10,10,2,0.1, 100, 10,
    20,25,4,0.3, 250, 30),
    byrow = TRUE,
    nrow = K)
  u_MRs <- matrix(c(
    10, 10, 2, 0.1, 100, 10,
    20, 25, 4, 0.3, 250, 30,  
    40, 50, 8, 0.5, 400, 200),
    byrow=TRUE,
    nrow = K)
  xi0 <- tau0 <- matrix(0, nrow = K, ncol = p)
  for (j in 1:p) {
    for (k in 1:K) {
      z <- interval_to_gaussian_general(l_MRs[k, j], u_MRs[k, j],0.99) # og=0.99
      xi0[k,j] <- z[1]
      tau0[k,j] <- z[2]
    } 
  }
  s0 <- matrix(0,nrow=L0,ncol=p)
  for (j in 1:p) {
    s0[,j] <- sample(3,3,T)
  }
  ### cluster centers and shapes ###
  mu0 <- matrix(0, nrow = L0, ncol = p) # true mean
  for (l in 1:L0) {
    for (j in 1:p) {
      mu0[l,j] <- rnorm(n=1,
                        mean = xi0[s0[l,j],j],
                        sd = sqrt(tau0[s0[l,j],j]))
    }
  }
  sigma0 <- matrix(c(1/6,1/2,0.6,0.05,20,20,
                     1/6,1/2,0.6,0.05,20,20,
                     1/6,1/2,1/3,0.05,20,20),
                   byrow=TRUE,
                   nrow = L0) # true variance
  x <- matrix(0, nrow = n, ncol = p)
  #nu0 <- 5
  for (i in 1:n) {
    for (j in 1:p) {
      x[i,j] <- mu0[c0[i],j] + sigma0[c0[i],j]*rnorm(n=1)
      ## extra options
      #rsn(n=1, xi = mu0[c0[i],j], omega=sigma0[c0[i],j], alpha=nu0[c0[i]]) # skew-normal
      #rnorm(n=1,  mean = mu0[c0[i],j], sd = sigma0[c0[i],j]) # normal
      #rnorm(n=1,  mean = mu0[c0[i],j], sd = sigma0[c0[i],j]) 
      #+ nu0*runif(1) # perturbed normal
    }
  }
  #pairs(x)
  ### fitting CLAMR
  # first adjust parameters
  pris <- write_sim_priors(p=p,K=K,l_MRs = l_MRs, u_MRs = u_MRs, xi0 = xi0, tau0 = tau0)
  mu.old.center = pris$mu.old.center
  mu.old.scale = pris$mu.old.scale
  var.priors = pris$var.priors
  source("rfuncts/marginal_cluster_list.R")
  source("rfuncts/tune_rho.R")
  R <- 15000
  Bu <- 5000
  L <- 10
  Kvec <- rep(K,p)
  stops <- R
  fit <- cluster.gibbs.topiclist(R=R,
                                 B=Bu,
                                 seed=seed,
                                 x=x,
                                 var.priors = var.priors,
                                 mu.old.center = mu.old.center,
                                 mu.old.scale = mu.old.scale,
                                 miss = is.na(x),
                                 K=Kvec,
                                 L=L,
                                 alpha_lambda = 1,
                                 rho = rep(0.7,p), 
                                 stops = stops)
  c.samps <- fit$c
  s.samps <- fit$s
  xhat <- fit$xhat
  # point estimate
  c.psm <- mcclust::comp.psm(c.samps)
  c.mv <- mcclust.ext::minVI(c.psm, c.samps)
  c.VI <- c.mv$cl
  # competitors
  ### Mclust EM Algorithm ###
  mcl <- Mclust(x, G=1:5)
  c.mcl <- mcl$classification
  ### Vanilla Mixture Model ###
  sc.x <- scale(x)
  fit.vanilla <- vanilla(R=R,
                         B=Bu,
                         seed=seed,
                         x=sc.x,
                         miss = is.na(x),
                         mu.priors = matrix(c(0,1), nrow = p, ncol = 2,byrow = TRUE),
                         var.priors = matrix(1, nrow = p, ncol = 2),
                         L=L,
                         gam = 1,
                         stops = stops)
  c.vanilla.samps <- fit.vanilla$c
  c.vanilla.psm <- mcclust::comp.psm(c.vanilla.samps)
  c.vanilla.mv <- mcclust.ext::minVI(c.vanilla.psm, c.vanilla.samps)
  c.vanilla <- c.vanilla.mv$cl
  ### k-means
  km <- kmeans(x,centers=L0)
  c.km <- km$cluster
  ### hiearachical clustering
  hca <- hclust(dist(x,method="euclidean"), method = "complete")
  c.hca <- cutree(hca,k=L0)
  ### Perfect Mixture Model ###
  fit.perfect <- vanilla(R=R,
                         B=Bu,
                         seed=seed,
                         x=x,
                         miss = is.na(x),
                         mu.priors = cbind(apply(mu0,2,mean), apply(mu0,2,var)),
                         var.priors = cbind(100,99*colMeans(sigma0^2)),
                         L=L0,
                         gam = 1,
                         stops = stops)
  c.perfect.samps <- fit.perfect$c
  c.perfect.psm <- mcclust::comp.psm(c.perfect.samps)
  c.perfect.mv <- mcclust.ext::minVI(c.perfect.psm, c.perfect.samps)
  c.perfect <- c.perfect.mv$cl
  # comparison
  ari <- c(
    adjustedRandIndex(c.VI,c0), # CLAMR
    adjustedRandIndex(c.vanilla,c0), # vanilla GMM
    adjustedRandIndex(c.mcl,c0), # mclust
    adjustedRandIndex(c.km,c0), # k-means
    adjustedRandIndex(c.hca,c0), # hca
    adjustedRandIndex(c.perfect,c0) # perfect
  )
  Kn <- c(    
    length(table(c.VI)), # CLAMR
    length(table(c.vanilla)), # vanilla GMM
    length(table(c.mcl)), # mclust
    length(table(c.km)), # k-means
    length(table(c.hca)), # hca
    length(table(c.perfect))# perfect)
  )
  return(list(ari=ari, Kn = Kn))
}