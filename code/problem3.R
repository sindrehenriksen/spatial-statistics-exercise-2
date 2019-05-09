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
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme(plot.title = element_text(hjust = 0.5))

# finding optimal number of clusters
nb <-fviz_nbclust(redwood_df, kmeans,
             method = "gap_stat")

redwood.opti.cluster.plot <- nb + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        axis.line = element_line())
redwood.opti.cluster.plot

ggsave("../figures/numb_clusters.pdf", plot = redwood.opti.cluster.plot, device = NULL, path = NULL,
       scale = 1, width = 5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)
# partitioning clusters k-medoids
km.res <- kmeans(redwood_df, 7, nstart = 25)
redwood.cluster.plot <-fviz_cluster(km.res, data = redwood_df,
                                    ellipse.type = "norm",
                                    palette = "Set1",
                                    stand = FALSE,
                                    geom = "point",
                                    ggtheme = theme_bw())
redwood.cluster.plot<- redwood.cluster.plot + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
# save plot
redwood.cluster.plot
ggsave("../figures/redwood_cluster_partitioning.pdf", plot = redwood.cluster.plot, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

redwood.cluster <- redwood.cluster.plot$data

# children position gaussian distribution
sigma_c <- mean(km.res$tot.withinss/km.res$size)

# number children poisson distributed
lambda_c <- mean(km.res$size)

# neuman-scott event RF
neuman_scott <- function(lambda_c,sigma_c,lambda_m){
  k <- 0 
  km = 0
  while (km == 0){
    km <- rpois(1,lambda_m)
  }
  x = numeric(500)
  y = numeric(500)
  cluster = numeric(500)
  for (j in seq(1,km)){
    tmp_x = runif(1)
    tmp_y = runif(1)
    kc = rbinom(1,length(redwood$x),1/km) #rpois(1,lambda_c)
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
set.seed(31)
x.ns<-neuman_scott(lambda_c,sigma_c,7)
cluster.event.plot <- ggplot(data = x.ns,aes(x=x, y=y)) + 
  geom_point(aes(colour = factor(cluster))) +
  scale_colour_manual(values=rep(cbPalette, length = length(x.ns$x)))+
  ggtitle("Clustered event RF") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(colour = "Mother")+
  xlim(0, 1)+
  ylim(0, 1)
cluster.event.plot
ggsave("../figures/cluster_event_rf.pdf", plot = cluster.event.plot, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


set.seed(15)
x.ns1 <- neuman_scott(lambda_c, sigma_c,7)
rel.cluster.plot1 <- ggplot(data = x.ns1,aes(x=x, y=y)) + 
  geom_point() +
  ggtitle("Realization") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())+
  xlim(0, 1)+
  ylim(0, 1)
rel.cluster.plot1

set.seed(5)
x.ns2 <- neuman_scott(lambda_c, sigma_c,7)
rel.cluster.plot2 <- ggplot(data = x.ns2,aes(x=x, y=y)) + 
  geom_point() +
  ggtitle("Realization") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())+
  xlim(0, 1)+
  ylim(0, 1)

rel.cluster.plot2

set.seed(7)
x.ns3 <- neuman_scott(lambda_c, sigma_c,7)
rel.cluster.plot3 <- ggplot(data = x.ns3,aes(x=x, y=y)) + 
  geom_point() +
  ggtitle("Realization") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())+
  xlim(0, 1)+
  ylim(0, 1)
rel.cluster.plot3

data.cluster.plot <- ggplot(data = redwood_df,aes(x=x,y=y)) + 
  geom_point() +
  ggtitle("Redwood data") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())+
  xlim(0, 1)+
  ylim(0, 1)
data.cluster.plot
rel.cluster.grid <- grid.arrange(
  rel.cluster.plot1,
  rel.cluster.plot2,
  rel.cluster.plot3,
  data.cluster.plot
)
ggsave("../figures/cluster_rel.pdf", plot = rel.cluster.grid, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------
# Monte Carlo test
set.seed(2)

#Condition on number of observed points

num.samps = 100
temp.rows = 1000
sample.ns.x <-matrix(NA,nrow = temp.rows,ncol = num.samps)
sample.ns.y <-matrix(NA,nrow = temp.rows,ncol = num.samps)
L.samps.ns.t = matrix(NA, nrow = 70, ncol = num.samps)
L.samps.ns.L = matrix(NA, nrow = 70, ncol = num.samps)

for (i in seq(1,num.samps)){
  temp.sample.cluster = neuman_scott(lambda_c,sigma_c,9)
  sample.ns.x[,i] = c(temp.sample.cluster$x,rep(NA,temp.rows-length(temp.sample.cluster$x)))
  sample.ns.y[,i] = c(temp.sample.cluster$y,rep(NA,temp.rows-length(temp.sample.cluster$y)))
  temp.sample.cluster.list = as.list(temp.sample.cluster[-3])
  temp.L.ns = Kfn(temp.sample.cluster.list, fs = 1)
  L.samps.ns.t[,i] = temp.L.ns$x
  L.samps.ns.L[,i] = temp.L.ns$y
}


ns.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  ns.quantiles[i,]= quantile(L.samps.ns.L[i,], c(0.05,0.95))
}

#Empirical 0.98-intervals
ns.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  ns.quantiles2[i,]= quantile(L.samps.ns.L[i,], c(0.01,0.99))
}

L_df = as.data.frame(L.samps.ns.L)
L_df$t = L.samps.ns.t[,1]
#Changing format to be able to plot all realisations
long_L = melt(L_df, id = 't')
L_redwood = Kfn(redwood, fs = 1) 
L_redwood_df = data.frame(t = L_redwood$x, L = L_redwood$y)
ns1.plot = ggplot(long_L,
                  aes(x=t, y=value, colour=variable)) +
  scale_colour_manual(values=rep(cbPalette, length = num.samps))+
  geom_line()+
  geom_line(data = ns.quantiles, aes(y = lower, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  geom_line(data = ns.quantiles, aes(y = upper, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  #geom_line(data = L_redwood_df, aes(y = L), col = 'black' )+
  ggtitle("Generated L-functions: Neuman-Scott") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

ns1.plot
ggsave("../figures/gen_ns_l.pdf", plot = ns1.plot, device = NULL, path = NULL,
       scale = 1, width = 6, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)

ns2.plot = ggplot(data = ns.quantiles) +
  geom_line(aes(y = lower, x = L_df$t), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t), col = 'black')+
  geom_line(data = L_redwood_df, aes(y = L, x = t), col = 'red' )+
  ggtitle("Redwood L-function with Neuman-Scott event RF 0.05 and 0.95 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

ns2.plot

ggsave("../figures/ns_quant1.pdf", plot = ns2.plot, device = NULL, path = NULL,
       scale = 1, width = 6, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)

t_max = 0.25
l = which(L_df$t ==t_max)
ns3.plot = ggplot(data = ns.quantiles2[1:l,]) +
  geom_line(aes(y = lower, x = L_df$t[1:l]), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t[1:l]), col = 'black')+
  geom_line(data = L_redwood_df[1:l,], aes(y = L, x = t), col = 'red' )+
  ggtitle("Redwood L-function with Neuman-Scott event RF 0.01 and 0.99 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

ns3.plot

ggsave("../figures/ns_quant2.pdf", plot = ns3.plot, device = NULL, path = NULL,
       scale = 1, width = 6, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)
