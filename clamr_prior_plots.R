library(ggplot2)
library(grid)
library(gridExtra)
#### cluster centers ####
set.seed(1996)
nsamp <- 100
x1 <- mvtnorm::rmvnorm(nsamp, mean = c(2,2), sigma = diag(c(1,1)))
x2 <- mvtnorm::rmvnorm(nsamp, mean = c(-2,-2), sigma = diag(c(1.5,1.5)))
x <- as.data.frame(rbind(x1,x2, c(2,2), c(-2,-2)))
x$cl <- as.factor(c(rep(1,nsamp), rep(2,nsamp), 3, 3))
x$size <- c(rep(1,2*nsamp), 2,2)
x$shape <- as.factor(c(rep(1,2*nsamp),3,3))
ggplot(x,aes(x=V1,y=V2,col=cl,size=size, shape=shape)) + geom_point() + theme_bw() +
  xlab(" ") + ylab(" ") +
  scale_size(range=c(2,5),
             guide="legend") +
  scale_color_brewer(palette="Set2") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size=15))

# CLAMR prior visualization
source("rfuncts/interval_to_gaussian.R")
range.mat <- matrix(c(-3,-1,
                      -1,0,
                      0,1.5), nrow = 3, ncol = 2,byrow=TRUE)
ov.range <- c(-3,1.5)
g.pars <- matrix(0, nrow = 3, ncol = 2)
for (i in 1:nrow(range.mat)) {
  g.pars[i,] <- interval_to_gaussian(range.mat[i,1], range.mat[i,2])
}
ov.pars <- interval_to_gaussian(ov.range[1],ov.range[2])
x <- seq(-3.5+mean(ov.range),3.5+mean(ov.range),by=0.01)
w <- c(1/3,1/3,1/3)
y1 <- w[1]*dnorm(x, mean = g.pars[1,1], sd = sqrt(g.pars[1,2])) +
  w[2]*dnorm(x,mean = g.pars[2,1], sd = sqrt(g.pars[2,2]) ) +
  w[3]*dnorm(x,mean = g.pars[3,1], sd = sqrt(g.pars[3,2]) )
y2 <- dnorm(x, mean = ov.pars[1], ov.pars[2])
y3 <- dnorm(x, mean = ov.pars[1], sqrt(4^2))
clamr.df <- data.frame(x=x,
                       y1=y1,
                       y2=y2,
                       y3=y3)
p1 <- ggplot(clamr.df, aes(x=x,y=y2)) + geom_point() +
  xlab("Variable") +
  ylab("Density") +
  geom_vline(xintercept = -3, color = "red") + 
  geom_vline(xintercept = -1, color = "red") + 
  geom_vline(xintercept = 0, color = "red") + 
  geom_vline(xintercept = 1.5, color = "red") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  labs(title = "Weakly Informative Prior")
p2 <- ggplot(clamr.df, aes(x=x,y=y1)) + geom_point() +
  xlab("Variable") +
  ylab("Density") +
  geom_vline(xintercept = -3, color = "red") + 
  geom_vline(xintercept = -1, color = "red") + 
  geom_vline(xintercept = 0, color = "red") + 
  geom_vline(xintercept = 1.5, color = "red") +
  theme(text = element_text(size=12))  +
  labs(title = "CLAMR Prior")
p3 <- ggplot(clamr.df, aes(x=x,y=y3)) + geom_point() +
  xlab("Variable") +
  ylab("Density") +
  geom_vline(xintercept = -3, color = "red") + 
  geom_vline(xintercept = -1, color = "red") + 
  geom_vline(xintercept = 0, color = "red") + 
  geom_vline(xintercept = 1.5, color = "red") +
  theme(text = element_text(size=12))  +
  labs(title = "Non-Informative")
grid.arrange(p3,p1,p2)
# all in one plot
clamr.oneplot <- data.frame(x = rep(x,3),
                            y = c(y1, y2, y3),
                            Prior = c(rep("CLAMR", length(x)), 
                                  rep("Weakly Informative", length(x)),
                                  rep("Non-Informative", length(x))))
ggplot(clamr.oneplot, aes(x=x,y=y,color=Prior)) + 
  geom_line(linewidth=1.5) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = -3, color = "purple", linetype="dashed") + 
  geom_vline(xintercept = -1, color = "purple", linetype="dashed") + 
  geom_vline(xintercept = 0, color = "purple", linetype="dashed") + 
  geom_vline(xintercept = 1.5, color = "purple", linetype="dashed") +
  xlab("Feature") +
  ylab("Density") +
  theme(text = element_text(size=12))  +
  labs(title = "Prior Choices") +
  theme_bw() +
  theme(text = element_text(size=13), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank())  
# black and white version
ggplot(clamr.oneplot, aes(x=x,y=y,group=Prior,linetype=Prior)) + 
  geom_line(linewidth=1.5) +
  scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = -3, linetype="dashed") + 
  geom_vline(xintercept = -1,  linetype="dashed") + 
  geom_vline(xintercept = 0,  linetype="dashed") + 
  geom_vline(xintercept = 1.5, linetype="dashed") +
  xlab("Feature") +
  ylab("Density") +
  theme(text = element_text(size=12))  +
  labs(title = "Prior Choices") +
  theme_bw() +
  theme(text = element_text(size=13), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank())  


