setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(factoextra)
library(fields)

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

phi <- function(tau,model.param){
  if ((tau>=0)&(tau<=model.param$tau_0)){
    phi = model.param$phi_0
  }else{
    phi = model.param$phi_0 * exp(-model.param$phi_1*(tau - model.param$tau_0))
  }
  return(phi)
}

acceptanceProb <- function(try_x,try_y,u,xd,yd,model.param){
  phi.sum = 0
  for (i in seq(1,length(xd))){
    phi.sum = phi.sum + 
              phi(sqrt((xd[u] - xd[i])^2 + (yd[u] - yd[i])^2), model.param) -
              phi(sqrt((try_x - xd[i])^2 + (try_y - yd[i])^2), model.param) 
  }
  return(min(1,exp(phi.sum)))
}

straussEventRF <- function(k,model.param){
  N = 150
  M = 300000
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  xd = runif(k)
  yd = runif(k)
  x0 = xd
  y0 = yd
  acceptance.prob = numeric(M)
  for (i in seq(1,M)){
    setTxtProgressBar(pb, i)
    u = round(runif(1,1,k))
    try_x = runif(1)
    try_y = runif(1)
    acceptance.prob[i] <- acceptanceProb(try_x,try_y,u,xd,yd,model.param)
    if (runif(1)<acceptance.prob[i]){
      xd[u] = try_x
      yd[u] = try_y
    }
  }
  return(list(df = data.frame(x=xd,y=yd,x0 = x0, y0 = y0),accep.prob = acceptance.prob))
}
k = length(cells_df$x)

model.param <- list(
  tau_0 = min(res.dist.fact),
  phi_0 = 70,
  phi_1 = 56
)
run_time <- system.time(repulsive.event.rf<-straussEventRF(k,model.param))
res.dist.fact.min <- min(get_dist(repulsive.event.rf$df[,1:2], stand = FALSE, method = "euclidean"))
ggplot(enframe(repulsive.event.rf$accep.prob))+
  geom_histogram(aes(x = value, y= ..density..),bins = 30)

ggplot(repulsive.event.rf$df) + 
  geom_point(aes(x=x,y=y), color ="red") 
  #geom_point(aes(x=x0,y=y0),color = "blue")


ggplot(enframe(repulsive.event.rf$accep.prob))+
  geom_point(aes(x = name, y= value))

