#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(NbClust)

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
#sigma_c <- mean(var(redwood_df$x[redwood.cluster$cluster==1]),
#                var(redwood_df$y[redwood.cluster$cluster==1]),
#                var(redwood_df$x[redwood.cluster$cluster==2]),
#                var(redwood_df$y[redwood.cluster$cluster==2]),
#                var(redwood_df$x[redwood.cluster$cluster==3]),
#                var(redwood_df$y[redwood.cluster$cluster==3]))
sigma_c <- km.res$tot.withinss/length(redwood_df$x)
# number children poisson distributed
temp_val <- as.data.frame(table(redwood.cluster$cluster))
lambda_c <- mean(temp_val$Freq)
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
  geom_point(aes(colour = factor(cluster)),size=2) +
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Clustered event RF") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(colour = "Cluster")+
  theme_minimal()+
  xlim(0, 1)+
  ylim(0, 1)
cluster.event.plot
ggsave("../figures/cluster_event_rf.pdf", plot = cluster.event.plot, device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
    
