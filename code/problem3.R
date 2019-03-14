#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

# read data
redwood = as.list(read.table('../data_files/redwood.dat', col.names = c('x', 'y')))

#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme(plot.title = element_text(hjust = 0.5))