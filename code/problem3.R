setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(reshape2)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(NbClust)
library(spatstat)

# read data
redwood = as.list(read.table('../data/redwood.dat', col.names = c('x', 'y')))

#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme(plot.title = element_text(hjust = 0.5))

# using 30 different method of determine number of clusters
# finding optimal number of clusters
pdf(file = NULL)
nb <- NbClust(redwood_df, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")
dev.off()
redwood.opti.cluster.plot <- fviz_nbclust(nb) + theme_minimal()
redwood.opti.cluster.plot
#save plot
ggsave("../figures/numb_clusters.pdf", plot = redwood.opti.cluster.plot, device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

# partitioning clusters k-medoids
km.res <- kmeans(redwood_df, 3, nstart = 3)
redwood.cluster.plot <-fviz_cluster(km.res, data = redwood_df,
                                    ellipse.type = "norm",
                                    palette = "Set1",
                                    stand = FALSE,
                                    geom = "point",
                                    ggtheme = theme_minimal())
redwood.cluster.plot
# save plot
ggsave("../figures/redwood_cluster_partitioning.pdf", plot = redwood.cluster.plot, device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

redwood.cluster <- redwood.cluster.plot$data

# children position gaussian distribution
sigma_c <- mean(km.res$withinss/(km.res$size-1))

# number children poisson distributed
lambda_c <- mean(km.res$size)

# neuman-scott event RF
neuman_scott <- function(lambda_c,sigma_c,lambda_m){
  k <- 0 
  km <- rpois(1,lambda_m)
  x = numeric(500)
  y = numeric(500)
  cluster = numeric(500)
  for (j in seq(1,km)){
    tmp_x = runif(1)
    tmp_y = runif(1)
    kc = rpois(1,lambda_c)
    for (i in seq(1,kc)){
      try_x = rnorm(1,tmp_x,sigma_c)
      try_y = rnorm(1,tmp_y,sigma_c)
      if ((try_x>=0)&(try_x<=1)&(try_y>=0)&(try_y<=1)){
        k = k + 1
        x[k] = try_x
        y[k] = try_y
        cluster[k] = j
      }
    }
  }
  return(data.frame(x = x[1:k],y = y[1:k],cluster=cluster[1:k]))
}

x.ns<-neuman_scott(lambda_c,sigma_c,3)
cluster.event.plot <- ggplot(data = x.ns,aes(x=x, y=y)) + 
  geom_point(aes(colour = factor(cluster))) +
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Clustered event RF") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(colour = "Mother")+
  theme_minimal()+
  xlim(0, 1)+
  ylim(0, 1)
cluster.event.plot
ggsave("../figures/cluster_event_rf.pdf", plot = cluster.event.plot, device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
    

L.cluster = Kfn(pp=ppp(x = x.ns$x, y=x.ns$y), fs = 1,k=100)
L.cluster.df = data.frame(t = L.cluster$x, L = L.cluster$y)

L.cluster.plot = ggplot(data = L.cluster.df, aes(x=t, y=L)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_path(col = 'red') + 
  ggtitle("L-interaction function for cluster effect") +
  theme(plot.title = element_text(hjust = 0.5))

L.cluster.plot

ggsave("../figures/L_cluster_plot.pdf", plot = L.cluster.plot, device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

# -----------------------------------------------------------------------------
# Monte Carlo test

set.seed(123)

#Condition on number of observed points
L.ns <-function(){
  num.samps = 100
  temp.rows = 200
  sample.ns.x <-matrix(NA,nrow = temp.rows,ncol = num.samps)
  sample.ns.y <-matrix(NA,nrow = temp.rows,ncol = num.samps)
  L.samps.ns.t = matrix(NA, nrow = temp.rows, ncol = num.samps)
  L.samps.ns.L = matrix(NA, nrow = temp.rows, ncol = num.samps)
  
  for (i in seq(1,num.samps)){
    temp.sample.cluster = neuman_scott(lambda_c,sigma_c,3)
    sample.ns.x[,i] = c(temp.sample.cluster$x,rep(NA,temp.rows-length(temp.sample.cluster$x)))
    sample.ns.y[,i] = c(temp.sample.cluster$y,rep(NA,temp.rows-length(temp.sample.cluster$y)))
    
    temp.L.ns = Kfn(pp = list(x = temp.sample.cluster$x, y = temp.sample.cluster$y), fs = 1)
    L.samps.ns.t = c(temp.L.ns$x,rep(NA,temp.rows-length(temp.L.ns$x)))
    L.samps.ns.L = c(temp.L.ns$x,rep(NA,temp.rows-length(temp.L.ns$y)))
  }
}

L.ns()
redwood.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  redwood.quantiles[i,]= quantile(L_samps_cells_L[i,], c(0.05,0.95))
}

#Empirical 0.98-intervals
redwood.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  redwood.quantiles2[i,]= quantile(L_samps_cells_L[i,], c(0.01,0.99))
}

L_df = as.data.frame(L_samps_redwood_L)
L_df$t = L_samps_redwood_t[,1]
#Changing format to be able to plot all realisations
long_L = melt(L_df, id = 't')

redwood1 = ggplot(long_L,
                  aes(x=t, y=value, colour=variable)) +
  scale_colour_manual(values=rep(cbPalette, length = Num_samps))+
  geom_line()+
  geom_line(data = redwood.quantiles, aes(y = lower, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  geom_line(data = redwood.quantiles, aes(y = upper, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  #geom_line(data = L_redwood_df, aes(y = L), col = 'black' )+
  ggtitle("Generated L-functions: Redwood") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

redwood2 = ggplot(data = redwood.quantiles) +
  geom_line(aes(y = lower, x = L_df$t), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t), col = 'black')+
  geom_line(data = L_redwood_df, aes(y = L, x = t), col = 'red' )+
  ggtitle("Redwood L-function with Poisson RF 0.05 and 0.95 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))