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
source("rfuncts/vanilla_gmm.R")
source("interval_prior_list.R")
source("rfuncts/simulate_clusters.R")
set.seed(1)
# parallel computing set-up
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
# note: original analysis requires INDITe datasets dat and clin (latter includes demographic vars)
# remove significant missingness
n_nas <- rowSums(is.na(dat))
dat <- dat[n_nas<7,]
clin <- clin[n_nas<7,]
# input missing data
mouse <- mice(dat)
dat.c <- complete(mouse)
# fitting
R <- 200000
Bu <- 10000 
stops <- R
L <- 10 
alpha_lambda <- 1
fit.vanilla <- vanilla(R=R,
                       B=Bu,
                       x=scale(dat.c),
                       miss = is.na(dat),
                       mu.priors = matrix(c(0,1), nrow = ncol(dat), ncol = 2,byrow = TRUE),
                       var.priors = matrix(1, nrow = ncol(dat), ncol = 2),
                       L=L,
                       gam = 1,
                       stops = stops,
                       seed=1)
c.vanilla.samps <- fit.vanilla$c
c.vanilla.psm <- mcclust::comp.psm(c.vanilla.samps)
c.vanilla.mv <- mcclust.ext::minVI(c.vanilla.psm, c.vanilla.samps)
c.vanilla <- c.vanilla.mv$cl
# plotting PSM
pairs(cbind(log(dat$urea), log(dat$crt), log(dat$totbili)),col=c.vanilla)
c.psm.ordered = c.vanilla.psm[order(as.numeric(as.character(c.vanilla))), 
                      rev(order(as.numeric(as.character(c.vanilla))))]
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
## demographics
table(c.vanilla) # general breakdown of clusters
# comparison with sex
round(diag(1/table(c.vanilla)) %*% table(c.vanilla, clin$sex),3)
# comparison with HIV status
round(diag(1/table(c.vanilla)) %*% table(c.vanilla, clin$hiv_class),3)
# comparison with advanced HIV
round(diag(1/table(c.vanilla)) %*% table(c.vanilla, clin$cd4pct<=14),3)
# comparison with malaria
round(diag(1/table(c.vanilla)) %*% table(c.vanilla, clin$malaria),3)
# comparison with outcome
round(diag(1/table(c.vanilla)) %*% table(c.vanilla, clin$inpatient_outcome),3)
# comparing to PE from CLAMR
#c.samps <- readRDS("indite_output/c_samps.rds")
c.psm <- mcclust::comp.psm(c.samps)
c.mv <- mcclust.ext::minVI(c.psm, c.samps)
c.VI <- c.mv$cl
adjustedRandIndex(c.VI, c.vanilla)
table(c.VI,c.vanilla)
c.VI.plot = dplyr::recode(as.factor(c.VI),
                          "2"="3",
                          "3"="2",
                          "4"="5",
                          "5"="4")
ggplot(data.frame(BUN = log(dat$urea), Creatinine = log(dat$crt), 
                  BGMM = as.factor(c.vanilla), CLAMR = c.VI.plot),
       aes(x=BUN, y=Creatinine,shape=CLAMR,col=BGMM)) + geom_point()
dat_x = prcomp(dat.c,scale=T)$x
ggplot(data.frame(PC1 = dat_x[,1], PC2 = dat_x[,2], 
                  BGMM = as.factor(c.vanilla), CLAMR = as.factor(c.VI)),
       aes(x=PC1, y=PC2,shape=BGMM,col=CLAMR,fill=CLAMR)) + geom_point(size=2) +
  scale_color_discrete(breaks=as.character(1:7)) + 
  scale_fill_discrete(breaks=as.character(1:7)) +
  scale_shape_manual(values=13:19) +
  theme_bw() +
  labs(title = "Comparison of CLAMR and BGMM") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=13))
### black and white version ###
bgmm_bw = ggplot(data.frame(PC1 = dat_x[,1], PC2 = dat_x[,2], 
                  BGMM = as.factor(c.vanilla), CLAMR = as.factor(c.VI)),
       aes(x=PC1, y=PC2,shape=BGMM)) + geom_point(size=2) +
  scale_color_discrete(breaks=as.character(1:7)) + 
  scale_fill_discrete(breaks=as.character(1:7)) +
  scale_shape_manual(values=13:19) +
  theme_bw() +
  labs(title = "BGMM Clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=13))
clamr_bw = ggplot(data.frame(PC1 = dat_x[,1], PC2 = dat_x[,2], 
                            BGMM = as.factor(c.vanilla), CLAMR = as.factor(c.VI)),
                 aes(x=PC1, y=PC2,shape=CLAMR)) + geom_point(size=2) +
  scale_color_discrete(breaks=as.character(1:7)) + 
  scale_fill_discrete(breaks=as.character(1:7)) +
  scale_shape_manual(values=13:19) +
  theme_bw() +
  labs(title = "CLAMR Clusters") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=13))
cowplot::plot_grid(clamr_bw,bgmm_bw,nrow=1,
                   labels=c("(a)","(b)"))
# boxplots
bps <- list()
for (j in 1:ncol(dat)) {
  bps[[j]] <- ggplot(data.frame(cluster = as.factor(c.vanilla),var = dat[,j]), aes(group=cluster, fill = cluster, y = var)) + 
    geom_boxplot() +
    theme_bw() + 
    theme(axis.text.y = element_blank()) +
    labs(title = colnames(dat)[j])
}
grid.arrange(bps[[1]], bps[[2]], bps[[3]], bps[[4]], bps[[5]],
             bps[[6]], bps[[7]], bps[[8]], bps[[9]], bps[[10]],
             bps[[11]], bps[[12]], bps[[13]], bps[[14]], bps[[15]],
             nrow = 3)
 