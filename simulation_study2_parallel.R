# libraries
library(mcclust)
library(mcclust.ext)
library(mclust)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
source("rfuncts/simulate_clusters_no_MRs.R")
# parallel computing set-up
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
## parameters
N <- c(100, 500, 750, 1000)
R <- 100
x <- foreach(i = 1:length(N)) %:%
  foreach(r = 1:R,
          .packages = c("sn", "mclust", "mcclust", "mcclust.ext"),
          .noexport = c("update_c_marginal")) %dopar% {
            simulate_clusters_no_MRs(n=N[i], seed=r)
          }
mat.list = list()
for (n in 1:length(N)) {
  mat.list[[n]]=matrix(unlist(x[[n]]), nrow = R, byrow = TRUE)
}
saveRDS(mat.list, "file name here")
# displaying
result <- matrix(c(
  colMeans(mat.list[[1]]),
  colMeans(mat.list[[2]]),
  colMeans(mat.list[[3]]),
  colMeans(mat.list[[4]])
), nrow = 4, byrow = T)
ari_result = result[,1:5]
L_result = result[,-(1:5)]
xtable::xtable(round(ari_result,3))
xtable::xtable(round(L_result,3))
sd.result <- matrix(c(
  apply(mat.list[[1]],2,sd),
  apply(mat.list[[2]],2,sd),
  apply(mat.list[[3]],2,sd),
  apply(mat.list[[4]],2,sd)
), nrow = 4, byrow = T)
