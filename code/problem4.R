#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


cell = as.list(read.table('../data/cells.dat', col.names = c('x', 'y')))