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
source("rfuncts/marginal_cluster_list.R")
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
# note: initial analysis requires INDITe datasets dat and clin (later includes demographic vars)
# remove significant missingness
n_nas <- rowSums(is.na(dat))
dat <- dat[n_nas<7,]
clin <- clin[n_nas<7,]
# input missing data
mouse <- mice(dat)
dat.c <- complete(mouse)

## method 1: EM algorithm
set.seed(1)
em <- Mclust(dat.c,G=1:10)
c.em <- em$classification

## method 2: k-means
### choose with elbow plot
set.seed(1)
max.k.km <- 10
tss <- rep(0,max.k.km)
for (k in 1:max.k.km) {
  tss[k] <- kmeans(dat.c,centers=k)$tot.withinss
}
plot(1:max.k.km,tss,main="K-Means Elbow Plot", xlab = "# of Clusters", ylab = "TWSS")
### choose k=4 clusters
km <- kmeans(dat.c,centers=4)
c.km <- km$cluster

## method 3: HCA
### first, we inspect dendrogram
set.seed(1)
dat.dist <- dist(dat.c,diag=TRUE)
hca <- hclust(dat.dist)
plot(hca)
### next, we will set the number of clusters equal to CLAMR
c.hca <- cutree(hca,k=6)

## method 4: CLAMR
c.clamr <- readRDS("indite_output/c_final.rds")

## now, combine the clusterings together
cluster.table <- data.frame(
  CLAMR = as.factor(c.clamr),
  EM = as.factor(c.em),
  kmeans = as.factor(c.km),
  HCA = as.factor(c.hca)
)
round(adjustedRandIndex(c.clamr,c.em),3)
round(adjustedRandIndex(c.clamr,c.km),3)
round(adjustedRandIndex(c.clamr,c.hca),3)



## plotting pcs
pca <- prcomp(scale(dat.c))
pca_x <- pca$x
cluster.table$PC1 <- pca_x[,1]
cluster.table$PC2 <- pca_x[,2]
clamr.plot <- ggplot(cluster.table,aes(x=PC1,y=PC2,color=CLAMR)) + geom_point() +
  theme_bw() +
  labs(title="(a) CLAMR",color="Cluster") +
  xlab("PC1") + ylab("PC2") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
em.plot <- ggplot(cluster.table,aes(x=PC1,y=PC2,color=EM)) + geom_point() +
  theme_bw() +
  labs(title="(b) EM",color="Cluster") +
  xlab("PC1") + ylab("PC2") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
kmeans.plot <- ggplot(cluster.table,aes(x=PC1,y=PC2,color=kmeans)) + geom_point() +
  xlab("PC1") + ylab("PC2") +
  labs(title="(c) K-Means",color="Cluster") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
hca.plot <- ggplot(cluster.table,aes(x=PC1,y=PC2,color=HCA)) + geom_point() +
  theme_bw() +
  labs(title="(d) HCA",color="Cluster") +
  xlab("PC1") + ylab("PC2") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
plot_grid(clamr.plot,em.plot,kmeans.plot,hca.plot,
          #labels=c("(a)","(b)","(c)","(d)"),
          nrow=2)

# black and white #
clamr.bw <- ggplot(cluster.table,aes(x=PC1,y=PC2,shape=CLAMR)) + geom_point() +
  theme_bw() +
  labs(title="(a) CLAMR",color="Cluster") +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
em.bw <- ggplot(cluster.table,aes(x=PC1,y=PC2,shape=EM)) + geom_point() +
  theme_bw() +
  labs(title="(b) EM",color="Cluster") +
  xlab("PC1") + ylab("PC2") + 
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
kmeans.bw <- ggplot(cluster.table,aes(x=PC1,y=PC2,shape=kmeans)) + geom_point() +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
  labs(title="(c) K-Means",color="Cluster") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
hca.bw <- ggplot(cluster.table,aes(x=PC1,y=PC2,shape=HCA)) + geom_point() +
  theme_bw() +
  labs(title="(d) HCA",color="Cluster") +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
plot_grid(clamr.bw,em.bw,kmeans.bw,hca.bw,
          #labels=c("(a)","(b)","(c)","(d)"),
          nrow=2)

