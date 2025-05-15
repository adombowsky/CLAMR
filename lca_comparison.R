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
library(poLCA)
library(GGally)
library(cowplot)
library(coda)
source("rfuncts/marginal_cluster_list.R")
source("interval_prior_list.R")
source("rfuncts/simulate_clusters.R")
set.seed(1)
# note: original analysis uses INDITe datasets dat and clin (latter has demographic vars)
# remove significant missingness
n_nas <- rowSums(is.na(dat))
dat <- dat[n_nas<7,]
clin <- clin[n_nas<7,]
# input missing data
mouse <- mice(dat)
dat.c <- complete(mouse)


### next, converting to categorical ###
ukey <- A.u[-c(4,8,14)]
lkey <- A.l[-c(4,8,14)]
n <- nrow(dat.c)
p <- ncol(dat.c)
dat.lca <- matrix(0,nrow=n,ncol=p)

for (i in 1:n) {
  for (j in 1:p) {
    if (dat.c[i,j] < lkey[[j]][1]) {
      next
    } else {
      for (k in 1:K[j]) {
        if (dat.c[i,j]<ukey[[j]][k]) {
          dat.lca[i,j] <- k
          break
        } else {
          dat.lca[i,j] <- K[j] + 1
        }
      }
    }
  }
}

### fitting the LCA model ###
dat.lca <- as.data.frame(dat.lca)
colnames(dat.lca) <- colnames(dat)
f<-cbind(crp,temperature,wbc,resp_rate,
         bicarb,pulse_rate,sbp,urea,crt,
         ast,totbili,plt,alb,na,gluc_rand) ~ 1
set.seed(1)
M <- poLCA(f, data = dat.lca,nclass=6,maxiter=10000)
c.M <- M$predclass

### comparing to CLAMR ###
#c.clamr <- readRDS("indite_output/c_final.rds")
### for even comparison, remove outlier from both ###
otl <- which(c.clamr==7)
c.clamr <- c.clamr[-otl]
c.M <- c.M[-otl]
adjustedRandIndex(c.clamr,c.M)
fossil::rand.index(c.clamr,c.M)
table(c.M,c.clamr)
### interpretations ###
prbs <- M$probs
## urea
prbs$urea # j=8
apply(prbs$urea,1,which.max)
prbs$crt # j=9
apply(prbs$crt,1,which.max)
prbs$ast # j=10
apply(prbs$ast,1,which.max)
prbs$totbili # j=11
apply(prbs$totbili,1,which.max)
