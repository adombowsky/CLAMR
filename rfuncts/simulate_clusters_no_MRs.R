simulate_clusters_no_MRs <- function(n, seed) {
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
    1,2,2,0.1, 100, 10,
    2,5,4,0.3, 250, 30),
    byrow = TRUE,
    nrow = K)
  u_MRs <- matrix(c(
    1, 2, 2, 0.1, 100, 10,
    2, 5, 4, 0.3, 250, 30,  
    4, 10, 8, 0.5, 400, 200),
    byrow=TRUE,
    nrow = K)
  xi0 <- tau0 <- matrix(0, nrow = K, ncol = p)
  for (j in 1:p) {
    for (k in 1:K) {
      z <- interval_to_gaussian_general(l_MRs[k, j], u_MRs[k, j],0.99)
      xi0[k,j] <- z[1]
      tau0[k,j] <- z[2]
    } 
  }
  ### simulate cluster means across the MRs (i.e. no MR structure)
  xi <- tau <- rep(0,p)
  for (j in 1:p) {
    z <- interval_to_gaussian_general(l_MRs[1,j], u_MRs[3,j],0.95)
    xi[j] <- z[1]
    tau[j] <- z[2]
  }
  ### cluster centers and shapes ###
  mu0 <- matrix(0, nrow = L0, ncol = p) # true mean
  for (l in 1:L0) {
    for (j in 1:p) {
      mu0[l,j] <- rnorm(n=1,
                        mean = xi[j],
                        sd = sqrt(tau[j]))
    }
  }
  sigma0 <- matrix(c(1/6,1/2,0.6,0.05,20,20,
                     1/6,1/2,0.6,0.05,20,20,
                     1/6,1/2,1/3,0.05,20,20),
                   byrow=TRUE,
                   nrow = L0) # true variance
  x <- matrix(0, nrow = n, ncol = p)
  nu0 <- 5
  for (i in 1:n) {
    for (j in 1:p) {
      x[i,j] <- mu0[c0[i],j] + sigma0[c0[i],j]*rt(n=1,nu0)
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
                         x=sc.x,
                         miss = is.na(x),
                         mu.priors = matrix(c(0,1), nrow = p, ncol = 2,byrow = TRUE),
                         var.priors = matrix(1, nrow = p, ncol = 2),
                         L=L,
                         gam = 1,
                         stops = stops,
                         seed=seed)
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
  # comparison
  ari <- c(
    adjustedRandIndex(c.VI,c0), # CLAMR
    adjustedRandIndex(c.vanilla,c0), # vanilla GMM
    adjustedRandIndex(c.mcl,c0), # mclust
    adjustedRandIndex(c.km,c0), # k-means
    adjustedRandIndex(c.hca,c0) # hca
  )
  Kn <- c(    
    length(table(c.VI)), # CLAMR
    length(table(c.vanilla)), # vanilla GMM
    length(table(c.mcl)), # mclust
    length(table(c.km)), # k-means
    length(table(c.hca)) # hca)
  )
  return(list(ari=ari, Kn = Kn))
}