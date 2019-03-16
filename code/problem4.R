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

findPhi <- function(distance){
  phi.res = numeric(length(distance))
  model.param = list(
    tau_0 = min(distance),
    phi_0 = 200,
    phi_1 = 60/(max(distance)-min(distance))
  )
  for (i in seq(1,length(distance))){
    phi.res[i] <- phi(distance[i],model.param)
  }
  phi.res = enframe(phi.res)
  phi.res.plot<- ggplot(phi.res,aes(x = name, y = value))+geom_point()
  return(list(phi.res = phi.res, phi.res.plot = phi.res.plot, model.param = model.param))
}

phi <- function(tau,model.param){
  if ((tau>=0)&(tau<model.param$tau_0)){
    phi = model.param$phi_0
  }else{
    phi = model.param$phi_0 * exp(-model.param$phi_1*(tau-model.param$tau_0))
  }
}
phi.list <- findPhi(res.dist.fact)

acceptanceProb <- function(try_x,try_y,u,xd,yd,model.param){
  phi.vec = numeric(length(xd))
  for (i in seq(1,length(xd))){
    phi.vec[i] = phi(sqrt((xd[u]-xd[i])^2+(yd[u]-yd[i])^2), model.param) - 
      phi(sqrt((try_x-xd[i])^2 + (try_y-yd[i])^2),model.param) 
  }
  return(min(1,exp(-sum(phi.vec))))
}

straussEventRF <- function(k,model.param){
  N = 150
  M = 300000
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  xd = runif(k)
  yd= runif(k)
  x0 = xd
  y0 = yd
  for (i in seq(1,M)){
    setTxtProgressBar(pb, i)
    u = sample(x=seq(1,k),size=1)
    try_x = runif(1)
    try_y = runif(1)
    acceptance.prob <- acceptanceProb(try_x,try_y,u,xd,yd,model.param)
    if (runif(1)<acceptance.prob){
      xd[u] = try_x
      yd[u] = try_y
    }
  }
  return(data.frame(x=xd,y=yd,x0 = x0, y0 = y0))
}
k = length(cells_df$x)

model.param <- list(
  tau_0 = median(res.dist.fact),
  phi_0 = 200,
  phi_1 = 90/(max(res.dist.fact)-min(res.dist.fact))
)
run_time <- system.time(repulsive.event.rf<-straussEventRF(k,model.param))
ggplot(repulsive.event.rf) + 
  geom_point(aes(x=x,y=y), color ="red") + 
  geom_point(aes(x=x0,y=y0),color = "blue")

