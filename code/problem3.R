#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(factoextra)
library(RColorBrewer)

# read data
redwood = as.list(read.table('../data/redwood.dat', col.names = c('x', 'y')))

#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme(plot.title = element_text(hjust = 0.5))

res.dist <- get_dist(redwood_df,stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_nbclust(redwood_df, kmeans, method = "gap_stat")
set.seed(123)
km.res <- kmeans(redwood_df, 7, nstart = 25)
a <-fviz_cluster(km.res, data = redwood_df,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())
