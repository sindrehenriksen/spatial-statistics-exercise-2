#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(factoextra)

cells = as.list(read.table('../data/cells.dat', col.names = c('x', 'y')))

cells_df = data.frame(x = cells$x, y = cells$y)

ggplot(data = cells_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Cells data") +
  theme(plot.title = element_text(hjust = 0.5))

res.dist <- rdist(cells_df)
res.dist.fact <- get_dist(cells_df, stand = FALSE, method = "euclidean")
fviz_dist(res.dist.fact, order = TRUE, show_labels = TRUE, lab_size = NULL,
          gradient = list(low = "red", mid = "white", high = "blue"))
findPhi <- function(res.dist.fact){
  phi_0 = 20
  tau_0 = min(res.dist)*0.6
  phi_1 = 6/(max(res.dist)-min(res.dist))
  phi.res = numeric(length(res.dist))
  for (i in seq(1,length(res.dist))){
    phi.res[i] <- phi(res.dist.vec[i],tau_0,phi_0,phi_1)
  }
  phi.res = enframe(phi.res)
  phi.res.plot<- ggplot(phi.test,aes(x = name, y = value))+geom_point()
  return(list(phi.res = phi.res, phi.res.plot = phi.res.plot, phi_0 = phi_0,phi_1 = phi_1, tau_0 = tau_0))
}

phi <- function(tau,tau_0,phi_0,phi_1){
  if ((tau>=0)&(tau<tau_0)){
    phi = phi_0
  }else{
    phi = phi_0 * exp(-phi_1*(tau-tau_0))
  }
}
phi.list <- findPhi(res.dist)
model.param <- list(tau_0 = phi.list$tau_0, phi_0 = phi.list$phi_0, phi_1 = phi.list$phi_1)
acceptanceProb <- function(try_x,try_y,u,xd,yd,model.param){
  phi.vec = numeric(size(xd,2))
  for (i in seq(1,size(xd,2))){
    
  }
  acceptance.prob = min(c(
    1,
    exp(-)
  ))
}

straussEventRF <- function(k,model.param){
  N = 10000
  xd = matrix(NA,nrow = N,ncol = k)
  yd = matrix(NA,nrow = N,ncol = k)
  xd[1,] = runif(k)
  yd[1,] = runif(k)
  for (i in seq(2,N)){
    xd[i,] = xd[i-1,]
    yd[i,] = yd[i-1,]
    u = runif(1,min = 0, max = 1)
    try_x = runif(1)
    try_y = runif(1)
    acceptance.prob <- acceptanceProb(try_x,try_y,u,xd,yd,model.param)
  }
}