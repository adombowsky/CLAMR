# libraries
library(ggplot2)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(mice)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(mcclust)
library(mcclust.ext)
library(mclust)
library(foreach)
library(doParallel)
library(GGally)
library(cowplot)
library(coda)
library(abind)
source("rfuncts/marginal_cluster_list.R")
source("rfuncts/simulate_clusters.R")
source("rfuncts/posterior_predictive.R")
set.seed(1)
### note: analysis requires INDITe datasets dat and clin (latter includes demographic vars)
# parallel computing set-up
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
# remove significant missingness
n_nas <- rowSums(is.na(dat))
dat <- dat[n_nas<7,]
clin <- clin[n_nas<7,]
# initialize missing data 
mouse <- mice(dat)
dat.c <- complete(mouse)
# fitting
R <- 300000 
B <- 10000 
stops <- R
L <- 10 
alpha_lambda <- 1
rho <- c()
n.chains <- 10
# setting rhoj
for (j in 1:ncol(dat)) {
  if (length(mu.old.center[[j]])==2) { # rho = 0.7 for Kj = 3, rho = 1.1 for Kj = 2
    rho[j] = 1.1
  } else {
    rho[j] = 0.7
  }
}
t1 <- Sys.time()
x <- foreach(r = 1:n.chains,
          .packages = c("sn", "mclust", "mcclust", "mcclust.ext"),
          .noexport = c("update_c_marginal")) %dopar% {
            cluster.gibbs.topiclist(R=R,
                                    B=B,
                                    seed=r,
                                    #x=as.matrix(dat.c),
                                    x = cbind(dat.c$urea,
                                              dat.c$crt,
                                              dat.c$ast,
                                              dat.c$totbili),
                                    #var.priors = var.priors,
                                     var.priors = var.priors[8:11,],
                                    #mu.old.center = mu.old.center,
                                     mu.old.center = mu.old.center[8:11],
                                    #mu.old.scale = mu.old.scale,
                                     mu.old.scale = mu.old.scale[8:11],
                                    #miss = is.na(dat),
                                    miss = is.na(cbind(dat$urea,
                                                 dat$crt,
                                                 dat$ast,
                                                 dat$totbili)),
                                    #K=K,
                                    K=K[8:11],
                                    L=L,
                                    alpha_lambda = alpha_lambda,
                                    #rho = rho,
                                     rho = rho[8:11],
                                    stops = stops)
          }
## analyzing the output
t2 <- Sys.time()
t.final <- round(as.numeric(difftime(t2,t1,units="secs")),3)
c.samps = c()
maxsize.samps = c()
#s.samps = c()
#ll.samps = c()
ind <- seq(1, nrow(x[[1]]$c), by = 10)
#pes <- matrix(0, nrow = n.chains, ncol = nrow(dat))
for (r in 1:n.chains) {
  # psm_r = mcclust::comp.psm(x[[r]]$c[ind,])
  # pes[r,] <- mcclust.ext::minVI(psm_r, x[[r]]$c[ind,])$cl
  maxsize.samps = cbind(maxsize.samps, apply(x[[r]]$c[ind,], 1, function(x) max(table(x))))
  # ll.samps = cbind(ll.samps, x[[r]]$ll[ind])
}
L.samps <- apply(c.samps,1,function(x) length(table(x)))
# maximum size traceplots
rstan::Rhat(maxsize.samps[,-5])
plot(maxsize.samps[,1],type="l", ylim = c(60,200))
lines(maxsize.samps[,2],col='red')
lines(maxsize.samps[,3],col='blue')
lines(maxsize.samps[,4],col='purple')
lines(maxsize.samps[,5],col='green')
lines(maxsize.samps[,6],col='grey')
lines(maxsize.samps[,7],col='pink')
lines(maxsize.samps[,8],col='brown')
lines(maxsize.samps[,9],col='orange')
lines(maxsize.samps[,10],col='cyan')
combinedchains = mcmc.list(as.mcmc(maxsize.samps[,1]), 
                           as.mcmc(maxsize.samps[,2]), as.mcmc(maxsize.samps[,3]), 
                           as.mcmc(maxsize.samps[,4]), 
                           #as.mcmc(maxsize.samps[,5]),
                           as.mcmc(maxsize.samps[,6]), as.mcmc(maxsize.samps[,7]), 
                           as.mcmc(maxsize.samps[,8]), as.mcmc(maxsize.samps[,9]), as.mcmc(maxsize.samps[,10]))
plot(combinedchains)
## compute overall ESS
maxsize.ov <- tidyr::gather(as.data.frame(maxsize.samps[,-5]))$value
coda::effectiveSize(maxsize.ov)
# log-likelihood traceplots
rstan::Rhat(ll.samps[,-5])
plot(ll.samps[,1],type="l", ylim = c(min(ll.samps), max(ll.samps)))
lines(ll.samps[,2],col='red')
lines(ll.samps[,3],col='blue')
lines(ll.samps[,4],col='purple')
lines(ll.samps[,5],col='green')
lines(ll.samps[,6],col='grey')
lines(ll.samps[,7],col='pink')
lines(ll.samps[,8],col='brown')
lines(ll.samps[,9],col='orange')
lines(ll.samps[,10],col='cyan')
apply(ll.samps,2,coda::effectiveSize)

combinedchains = mcmc.list(as.mcmc(ll.samps[,1]), as.mcmc(ll.samps[,2]), as.mcmc(ll.samps[,3]), 
                           as.mcmc(ll.samps[,4]), 
                           as.mcmc(ll.samps[,5]),
                           as.mcmc(ll.samps[,6]), as.mcmc(ll.samps[,7]), as.mcmc(ll.samps[,8]), as.mcmc(ll.samps[,9]), as.mcmc(ll.samps[,10]))
plot(combinedchains)
ll.ov <- tidyr::gather(as.data.frame(ll.samps[,-5]))$value
coda::effectiveSize(ll.ov)
# comparing point estimates
pes <- pes[-5,]
D <- matrix(0, nrow = nrow(pes), ncol = nrow(pes))
for (i in 1:nrow(pes)) {
  for (j in 1:nrow(pes)) {
    D[i,j] <- fossil::rand.index(pes[i,], pes[j,])
  }
}
summary(D[upper.tri(D)])
# using ggplot
melted_df_ll = reshape2::melt(ll.samps[,-5], id.vars = 2)
LL_chains = ggplot(data.frame(Chain = as.factor(melted_df_ll$Var2), LL = melted_df_ll$value, Iteration = 1:nrow(ll.samps)), 
       aes(x = Iteration, y = LL, color = Chain)) +
  geom_line(alpha=0.8) +
  labs(title = "Log-Likelihood") +
  theme_bw() +
  theme(title = element_text(size = 12))
melted_df_max = reshape2::melt(maxsize.samps[,-5], id.vars = 2)
max_chains = ggplot(data.frame(Chain = as.factor(melted_df_max$Var2), nmax = melted_df_max$value, Iteration = 1:nrow(maxsize.samps)), 
                   aes(x = Iteration, y = nmax, color = Chain)) +
  geom_line(alpha=0.8) +
  labs(title = latex2exp::TeX("$n^{\\max}$")) +
  theme_bw() +
  theme(title = element_text(size = 12))
leg <- get_legend(LL_chains)
plot_grid(LL_chains + theme(legend.position = "none"), 
          max_chains + theme(legend.position = "none"),
          leg,
          nrow = 1,
          rel_widths = c(1,1,0.2))
melted_df_rand = reshape2::melt(rand.samps[,-5], id.vars = 2)
rand_chains = ggplot(data.frame(Chain = as.factor(melted_df_rand$Var2), nmax = melted_df_rand$value, Iteration = 1:nrow(rand.samps)),
                    aes(x = Iteration, y = nmax, color = Chain)) +
  geom_line(alpha=0.8) +
  labs(title = latex2exp::TeX("$R(\\textbf{c}, \\textbf{c}^*)$")) +
  theme_bw() +
  theme(title = element_text(size = 12))
leg <- get_legend(rand_chains)
plot_grid(max_chains + theme(legend.position = "none"),
          rand_chains + theme(legend.position="none"),
          leg,
          nrow = 1,
          rel_widths = c(1,1,0.2))

# testing feature influence
# R.b <- dim(s.samps)[1]
# b <- bf <- c()
# epsilon <- 0.1
# for (j in 1:ncol(dat)) {
#   for (r in 1:R.b) {
#     levs <- as.numeric(levels(as.factor(c.samps[r,]))) # filled clusters
#     s.r <- s.samps[r, levs, j]
#     b[r] <- Ln_binder(s.r)
#   }
#   bf[j] <- (1-mean(b<epsilon))/mean(b<epsilon)
# }
# names(bf) <- colnames(dat)
# round(bf,3)

# computing the point estimate, skip chain # 4
c.samps = mu.samps = sigma_sq.samps =  c()
for (r in 1:n.chains) {
  if (r==5) { ## og: r=4
    next
  } else {
    c.samps = rbind(c.samps,x[[r]]$c[ind,])
    mu.samps = abind(mu.samps,x[[r]]$mu,along=3)
    sigma_sq.samps = abind(sigma_sq.samps,x[[r]]$sigma_sq,along=3)
  }
}
c.psm <- mcclust::comp.psm(c.samps)
c.mv <- mcclust.ext::minVI(c.psm, c.samps)
c.VI <- c.mv$cl
# recoding by cluster size
c.VI <- as.numeric(as.character((dplyr::recode(as.factor(c.VI),
                          "2"="4",
                          "3"="2",
                          "4"="3"))))
## demographics
table(c.VI) # general breakdown of clusters
# comparison with sex
round(diag(1/table(c.VI)) %*% table(c.VI, clin$sex),3)
# comparison with HIV status
round(diag(1/table(c.VI)) %*% table(c.VI, clin$hiv_class),3)
# comparison with advanced HIV
round(diag(1/table(c.VI)) %*% table(c.VI, clin$cd4pct<=14),3)
# comparison with malaria
round(diag(1/table(c.VI)) %*% table(c.VI, clin$malaria_rtest),3)
# comparison with outcome
round(diag(1/table(c.VI)) %*% table(c.VI, clin$inpatient_outcome),3)
# age
round(mean(clin$age),1)
round(sd(clin$age),1)
round(mean(clin$age[c.VI==1]),1)
round(sd(clin$age[c.VI==1]),1)
round(mean(clin$age[c.VI==2]),1)
round(sd(clin$age[c.VI==2]),1)
round(mean(clin$age[c.VI==3]),1)
round(sd(clin$age[c.VI==3]),1)
round(mean(clin$age[c.VI==4]),1)
round(sd(clin$age[c.VI==4]),1)
round(mean(clin$age[c.VI==5]),1)
round(sd(clin$age[c.VI==5]),1)
round(mean(clin$age[c.VI==6]),1)
round(sd(clin$age[c.VI==6]),1)

# plotting
# identify outlier
otl <- which(c.VI==7)
c.VI.plot = as.factor(c.VI)
x_r <- data.frame(BUN = log(dat$urea),
             Creatinine = log(dat$crt),
             AST = log(dat$ast),
             Bilirubin = log(dat$totbili),
             Cluster = c.VI.plot)
GGally::ggpairs(x_r[-otl,], columns = 1:4, aes(color=as.factor(Cluster)), upper="blank", legend=1) + 
  labs(fill="Cluster", title="Clusters and Influential Features") +
  theme_bw() +
  scale_fill_discrete(breaks=as.character(1:7)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=13))
### black and white version ###
## to get legend
ggplot(x_r[-otl,], aes(x=BUN, linetype = as.factor(Cluster))) + geom_density() + 
  labs(linetype = "Cluster")
ggplot(x_r[-otl,], aes(x=BUN, y=Creatinine, shape = as.factor(Cluster))) + geom_point() + 
  labs(shape = "Cluster")
GGally::ggpairs(x_r[-otl,], columns = 1:4, mapping=aes(group=as.factor(Cluster), 
                                               linetype=as.factor(Cluster),
                                               shape = as.factor(Cluster)
                                               ), 
                upper="blank", legend=1) + 
  labs(linetype="Cluster", 
       shape = "Cluster",
       title="Clusters and Influential Features") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=13))
plots <- list()
x_plot <- exp(x_r[,-5])
x_plot <- x_plot %>%
  mutate(BUN = BUN*0.357,
       Creatinine = Creatinine*88.4,
       Bilirubin =Bilirubin*17.104)#%>%
K.n <- length(table(c.VI))
for (j in 1:ncol(x_plot)) {
  plot.df <- data.frame(var = x_plot[,j], cluster = as.integer(c.VI.plot))
  plots[[j]] <- ggplot(plot.df[-otl,], aes(x=cluster, y=var, group = as.factor(cluster),  fill = as.factor(cluster))) + 
    geom_boxplot() + 
    scale_fill_discrete(breaks=as.character(1:6)) +
    labs(title = colnames(x_plot)[j], fill = "Cluster") +
    ylab(" ") + xlab(" ") +
    theme_bw()  
}
plot_grid(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], 
          nrow = 2)
# boxplots of characteristics outside of influential features
bps <- list()
box_titles <- c("C-Reactive Protein", "Temperature", "WBC", "Respiratory Rate",
                "Bicarbonate", "Pulse Rate", "SBP", "BUN", "Creatinine",
                "AST", "Bilirubin", "Platelets", "Albumin", "Sodium",
                "Glucose")
for (j in 1:ncol(dat)) {
  bps[[j]] <- ggplot(data.frame(cluster = as.integer(c.VI.plot),var = dat[,j])[-otl,],
                     aes(group=as.factor(cluster), fill = as.factor(cluster), 
                         x = cluster, y = var)) + 
    geom_boxplot() +
    theme_bw() + 
    ylab(" ") +
    xlab(" ") +
    theme(axis.text.y = element_blank()) +
    labs(title = box_titles[j], fill = "Cluster")
}
leg = get_legend(bps[[1]] +theme(legend.position="bottom"))
plot_grid(bps[[1]], 
             bps[[2]], 
             bps[[3]], 
             bps[[4]], 
             bps[[5]],
             bps[[6]],
             bps[[7]], 
             bps[[8]],
             bps[[9]],
             bps[[10]],
             bps[[11]], 
             bps[[12]], 
             bps[[13]], 
             bps[[14]], 
             bps[[15]],
             nrow = 3, ncol = 5)
# estimating ELs
x_plot_otl <- x_plot[-otl,]
c.VI.otl <- c.VI[-otl]
## AST
k <- 5
summary(x_plot_otl$AST[c.VI==k])
summary(x_plot_otl$Bilirubin[c.VI==k])
summary(x_plot_otl$BUN[c.VI==k])
summary(x_plot_otl$Creatinine[c.VI==k])
# credible ball
# credible ball
cb <- mcclust.ext::credibleball(c.VI, c.samps, c.dist = c("VI"))
c.horiz <- cb$c.horiz[1,]
x_cb <- data.frame(BUN = log(dat$urea),
                  Creatinine = log(dat$crt),
                  AST = log(dat$ast),
                  Bilirubin = log(dat$totbili),
                  Cluster = as.factor(c.horiz))
GGally::ggpairs(x_cb, columns = 1:4, aes(color=as.factor(Cluster)), upper="blank", legend=1) + 
  labs(fill="Cluster", title="Clusters and Influential Features") +
  theme_bw() +
  scale_fill_discrete(breaks=as.character(1:length(table(c.horiz)))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=13))
# centroid finding
x_red <- cbind(dat.c$urea,
               dat.c$crt,
               dat.c$ast,
               dat.c$totbili)
distvec = dist(x_red)
distmat = matrix(0, nrow = nrow(x_red), ncol = nrow(x_red))
distmat[upper.tri(distmat)] = distvec
distmat = distmat + t(distmat)
meds <- matrix(0, nrow = length(table(c.VI)), ncol = ncol(x_red))
for (j in 1:ncol(x_red)) {
  distvec = dist(x_red[,j])
  distmat = matrix(0, nrow = nrow(x_red), ncol = nrow(x_red))
  distmat[upper.tri(distmat)] = distvec
  distmat = distmat + t(distmat)
  for (l in 1:length(table(c.VI))) {
    meds[l,j] <- monoClust::medoid(which(c.VI==l), distmat)
  }
}
meds = meds[order(table(c.VI), decreasing = T),]
uq = medoid_prob(meds = meds, c_samps = c.samps, c = c.VI.plot)
ggplot(data.frame(v1 = log(x_red[,3]), v2 = log(x_red[,1]), prob = uq[,3]), aes(x=v1, y = v2, col = prob)) + geom_point(size=3)
# plotting PSM
c.psm.ordered = c.psm[order(as.numeric(as.character(c.VI.plot[-otl]))), 
                      rev(order(as.numeric(as.character(c.VI.plot[-otl]))))]
psm.df <- reshape2::melt(c.psm.ordered, c("x", "y"), value.name = "Pr")
ggplot(data=psm.df,aes(x=x,y=y,fill=Pr)) +
  geom_tile() +
  theme_bw() +
  scale_fill_distiller(palette = "RdPu") +
  labs(fill = latex2exp::TeX("$\\Pr(c_i = c_j | \\textbf{y})$"),
       title = "Posterior Similiarity Matrix") +
  xlab(" ") + ylab(" ") +
  theme(axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
### black and white version ### 
ggplot(data=psm.df,aes(x=x,y=y,fill=Pr)) +
  geom_tile() +
  theme_bw() +
  scale_fill_distiller(palette = "RdPu") +
  labs(fill = latex2exp::TeX("$\\Pr(c_i = c_j | \\textbf{y})$"),
       title = "Posterior Similiarity Matrix") +
  xlab(" ") + ylab(" ") +
  scale_fill_gradient(high = "white", low = "black") +
  theme(axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
# constructing WAIC
p_star <- l_squared <- l_dens <- c()
for (r in 1:n.chains) {
  if (r==5) { ## og: r=4
    next
  } else {
  p_star <- rbind(p_star,exp(x[[r]]$ldens[ind,]))
  l_squared <- rbind(l_squared,x[[r]]$ldens[ind,]^2)
  l_dens <- rbind(l_dens,x[[r]]$ldens[ind,])
  }
}
B_n = -(1/nrow(dat)) * sum(log(colMeans(p_star)))
V_n <- (1/nrow(dat))*(sum(colMeans(l_squared)) - sum(colMeans(l_dens)^2))
waic = B_n + V_n

## posterior predictive checks
x_pp <- posterior_predictive(c.samps,mu.samps,sigma_sq.samps,stops=10000)
x_pp_hat <- apply(x_pp$x, c(2,3), mean)
x_pp_low <- apply(x_pp$x,c(2,3),quantile,prob=0.025)
x_pp_high <- apply(x_pp$x,c(2,3),quantile,prob=0.975)

# subset equally spaced chains
x_pp_subset <- x_pp$x[seq(1,261000,by=29000),,]
urea_pp_subset <- x_pp_subset[,,1]
chains_df <- data.frame(urea = c(dat.c$urea,as.vector(t(urea_pp_subset))),
                        type = c(rep("Raw Data",nrow(dat.c)),
                                 rep("Chain 1", nrow(dat.c)),
                                 rep("Chain 2", nrow(dat.c)),
                                 rep("Chain 3", nrow(dat.c)),
                                 rep("Chain 4", nrow(dat.c)),
                                 rep("Chain 5", nrow(dat.c)),
                                 rep("Chain 6", nrow(dat.c)),
                                 rep("Chain 7", nrow(dat.c)),
                                 rep("Chain 8", nrow(dat.c)),
                                 rep("Chain 9", nrow(dat.c))
                                 ))
ggplot(chains_df,aes(x=urea,group=type,color=type)) + geom_density()

# overlay a bunch for each chain
its <- c(seq(1,nrow(x_pp$x),by=nrow(x_pp$x)/9), nrow(x_pp$x)+1)
r = 2 # chain
x_pp_subset <- x_pp[its[r]:(its[r+1]-1),,]
# generate 10 samples
n_samps <- 200
j <- 3 # variable
var_names <- c("BUN", "Creatinine", "AST", "Bilirubin")
smaller_subset <- x_pp_subset[seq(1,nrow(x_pp_subset),by=nrow(x_pp_subset)/n_samps),,]
#smaller_subset <- x_pp_subset[10:100,,]
chains_df <- data.frame(var = c(as.vector(t(smaller_subset[,,j])),x_plot[,j]),
                        type = c(rep(1:n_samps,each=nrow(dat.c)),
                                 rep("Raw Data",nrow(dat.c))),
                        name = c(rep("Samples",n_samps*nrow(dat.c)),
                                 rep("Raw Data",nrow(dat.c))),
                        width = c(rep(.1,n_samps*nrow(dat.c)),
                                  rep(.1,nrow(dat.c)))
                        )
ggplot(chains_df,aes(x=var,group=type,color=name)) + geom_density() +
  scale_color_manual(values=c("#EB4039","#D0BAEB")) +
  labs(title = paste("Chain",r, ",", var_names[j] )) +
  theme_bw()

# overlay a bunch for each chain
# generate 10 samples
n_samps <- 500
j <- 4 # variable
var_names <- c("BUN", "Creatinine", "AST", "Bilirubin")
smaller_subset <- x_pp[seq(1,nrow(x_pp),by=nrow(x_pp)/n_samps),,]
pp_plots <- list()
for (j in 1:4) {
  chains_df <- data.frame(var = c(as.vector(t(smaller_subset[,,j])),x_plot[,j]),
                          type = c(rep(1:n_samps,each=nrow(dat.c)),
                                   rep("Raw Data",nrow(dat.c))),
                          name = c(rep("Samples",n_samps*nrow(dat.c)),
                                   rep("Raw Data",nrow(dat.c))),
                          width = c(rep(.1,n_samps*nrow(dat.c)),
                                    rep(.1,nrow(dat.c)))
  )
  pp_plots[[j]] = ggplot(chains_df,aes(x=var,group=type,color=name)) + geom_density() +
    scale_color_manual(values=c("#EB4039","#D0BAEB")) +
    labs(title = paste(var_names[j] )) +
    xlab("Feature") +
    ylab("Density") +
    theme_bw() +
    theme(legend.position="none",
          text=element_text(size=12))
}
cowplot::plot_grid(pp_plots[[1]],pp_plots[[2]],pp_plots[[3]],pp_plots[[4]],
                   nrow=2, ncol=2)

# addressing posterior uncertainty in MR assignment
s.samps = c()
for (r in 1:n.chains) {
  if (r==1) {
    s.samps = x[[r]]$s[ind,,]
  } 
  else if (r==5) {
    next
  }
  else {
    s.samps = abind::abind(s.samps,x[[r]]$s[ind,,],along=1)
  }
}
prb_MR <- list()
K_inf <- K[8:11]
n_clusters <- length(table(c.VI))
K_j_max <- 3
# initialization
for (l in 1:n_clusters) {
  prb_MR[[l]] = array(0,dim=c(nrow(s.samps),K_j_max,4))
}
for (l in 1:n_clusters) {
  for (j in 1:4) {
    for (v in 1:K_inf[j]) {
      for (h in 1:nrow(s.samps)) {
        prb_MR[[l]][h,v,j] = mean(s.samps[h,c.samps[h,c.VI==l],j]==v)
      } 
    }
  } 
}
# point estimation
avg_prb_MR <- array(0,dim=c(n_clusters,4,K_j_max))
for (l in 1:n_clusters) {
  for (j in 1:4) {
    avg_prb_MR[l,j,] <- colMeans(prb_MR[[l]][,,j]) # clusters x variable x MRs
  }
}
max_prb_MR <- which_max_prb_MR <-  matrix(0,nrow=n_clusters,ncol=4)
for (l in 1:n_clusters) {
  for (j in 1:4) {
    max_prb_MR[l,j] <- max(avg_prb_MR[l,j,])
    which_max_prb_MR[l,j] <- which.max(avg_prb_MR[l,j,])
  }
}
# plotting posterior distribution
prbs_plots <- list()
prbs_names <- c("BUN", "CRT", "AST", "BRN")
for (l in 1:n_clusters) {
  prbs_plots[[l]] <- list()
  for (j in 1:4) {
    prbs_plots[[l]][[j]] <- ggplot(
      data = data.frame(x=prb_MR[[l]][,which_max_prb_MR[l,j],j]),
      aes(x=x)) + geom_histogram(color="blue",fill="skyblue1") +
      labs(title = paste("Cl", l, ",", prbs_names[j], ", MR", which_max_prb_MR[l,j])) +
      theme_bw() +
      xlab("Empirical Prob.") +
      theme(axis.text.y=element_blank(),
            axis.text.x = element_text(size=6.5),
            title = element_text(size=7))
  }
}

plot_grid(prbs_plots[[1]][[1]], prbs_plots[[1]][[2]], prbs_plots[[1]][[3]], prbs_plots[[1]][[4]],
          prbs_plots[[2]][[1]], prbs_plots[[2]][[2]], prbs_plots[[2]][[3]], prbs_plots[[2]][[4]],
          prbs_plots[[3]][[1]], prbs_plots[[3]][[2]], prbs_plots[[3]][[3]], prbs_plots[[3]][[4]],
          prbs_plots[[4]][[1]], prbs_plots[[4]][[2]], prbs_plots[[4]][[3]], prbs_plots[[4]][[4]],
          prbs_plots[[5]][[1]], prbs_plots[[5]][[2]], prbs_plots[[5]][[3]], prbs_plots[[5]][[4]],
          prbs_plots[[6]][[1]], prbs_plots[[6]][[2]], prbs_plots[[6]][[3]], prbs_plots[[6]][[4]],
          nrow=6
          )
