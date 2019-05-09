setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(spatial)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(factoextra)
library(fields)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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
  if (tau<=model.param$tau_0){
    phi.res = model.param$phi_0
  }else{
    phi.res = model.param$phi_0 * exp(-model.param$phi_1*(tau - model.param$tau_0))
  }
  return(phi.res)
}

acceptanceProb <- function(try_x,try_y,u,xd,yd,model.param){
  prev.sum = 0
  new.sum = 0
  for (i in seq(1,length(xd))){
    if (i!=u){
      prev.sum = prev.sum + phi(sqrt((xd[u] - xd[i])^2 + (yd[u] - yd[i])^2), model.param)
      new.sum = new.sum + phi(sqrt((try_x - xd[i])^2 + (try_y - yd[i])^2), model.param) 
    }
  }
  phi.sum = prev.sum-new.sum
  return(min(1,exp(phi.sum)))
}
set.seed(123)
straussEventRF <- function(k,model.param){
  M = 1000
  xd = runif(k)
  yd = runif(k)
  x0 = xd
  y0 = yd
  x1 = c(xd[1])
  y1 = c(yd[1])
  trace.dist = c(min(dist(data.frame(x=xd,y=yd))))
  acceptance.prob = numeric(M)
  for (i in seq(1,M)){
    u = round(runif(1,1,k))
    try_x = runif(1)
    try_y = runif(1)
    acceptance.prob[i] = acceptanceProb(try_x,try_y,u,xd,yd,model.param)
    if (runif(1)<acceptance.prob[i]){
      xd[u] = try_x
      yd[u] = try_y
      trace.dist = c(trace.dist, min(dist(data.frame(x=xd,y=yd))))
      x1 = c(x1,xd[1])
      y1 = c(y1,yd[1])
    }
  }
  return(list(df = data.frame(x=xd,y=yd,x0 = x0, y0 = y0),accep.prob = acceptance.prob,trace.dist = trace.dist,df2 = data.frame(x = x1,y = y1)))
}
k = length(cells_df$x)

model.param <- list(
  tau_0 = min(res.dist.fact),
  phi_0 = 7,
  phi_1 = 90
)
repulsive.event.rf<-straussEventRF(k,model.param)
res.dist.fact.min <- min(get_dist(repulsive.event.rf$df[,1:2], stand = FALSE, method = "euclidean"))
ggplot(enframe(repulsive.event.rf$accep.prob))+
  geom_histogram(aes(x = value, y= ..density..),bins = 30)

repulsive.event.plot<-ggplot(repulsive.event.rf$df) + 
  geom_point(aes(x=x,y=y))+
  ggtitle("Strauss Event RF") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
repulsive.event.plot
ggsave("../figures/repulsive_event_rf.pdf", plot = repulsive.event.plot, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

ggplot(enframe(repulsive.event.rf$accep.prob))+
  geom_point(aes(x = name, y= value))

ggplot(repulsive.event.rf$df2,aes(x=x,y=y)) + 
  geom_line()

repulsive.trace.plot<-ggplot(enframe(repulsive.event.rf$trace.dist)) + 
  geom_line(aes(x= name, y = value)) + 
  geom_hline(yintercept = model.param$tau_0, color = "red")+
  ggtitle("Trace plot of minimum distance") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        axis.line = element_line(),axis.title.x = element_blank())+ 
  scale_y_continuous(name = "minimum distance")
repulsive.trace.plot
ggsave("../figures/repulsive_trace_rf.pdf", plot = repulsive.trace.plot, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

ggplot(repulsive.event.rf$df) + 
  geom_point(aes(x=x0,y=y0))
set.seed(2)

#Condition on number of observed points

num.samps = 100
temp.rows = k
sample.s.x <-matrix(NA,nrow = temp.rows,ncol = num.samps)
sample.s.y <-matrix(NA,nrow = temp.rows,ncol = num.samps)
L.samps.s.t = matrix(NA, nrow = 70, ncol = num.samps)
L.samps.s.L = matrix(NA, nrow = 70, ncol = num.samps)

for (i in seq(1,num.samps)){
  temp.sample.repulsive = straussEventRF(k,model.param)
  sample.s.x[,i] = temp.sample.repulsive$df$x
  sample.s.y[,i] = temp.sample.repulsive$df$y
  temp.sample.repulsive.list = as.list(temp.sample.repulsive$df[1:2])
  temp.L.s = Kfn(pp = temp.sample.repulsive.list, fs = 1)
  L.samps.s.t[,i] = temp.L.s$x
  L.samps.s.L[,i] = temp.L.s$y
}


s.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  s.quantiles[i,]= quantile(L.samps.s.L[i,], c(0.05,0.95))
}

#Empirical 0.98-intervals
s.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  s.quantiles2[i,]= quantile(L.samps.s.L[i,], c(0.01,0.99))
}

L_df = as.data.frame(L.samps.s.L)
L_df$t = L.samps.s.t[,1]
#Changing format to be able to plot all realisations
long_L = melt(L_df, id = 't')
L_cells = Kfn(cells, fs = 1) 
L_cells_df = data.frame(t = L_cells$x, L = L_cells$y)
s1.plot = ggplot(long_L,
                  aes(x=t, y=value, colour=variable)) +
  scale_colour_manual(values=rep(cbPalette, length = num.samps))+
  geom_line()+
  geom_line(data = s.quantiles, aes(y = lower, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  geom_line(data = s.quantiles, aes(y = upper, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  #geom_line(data = L_redwood_df, aes(y = L), col = 'black' )+
  ggtitle("Generated L-functions: Strauss") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

s1.plot
ggsave("../figures/gen_strauss_l.pdf", plot = s1.plot, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)

s2.plot = ggplot(data = s.quantiles) +
  geom_line(aes(y = lower, x = L_df$t), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t), col = 'black')+
  geom_line(data = L_cells_df, aes(y = L, x = t), col = 'red' )+
  ggtitle("Cells L-function with Strauss event RF 0.05 and 0.95 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

s2.plot

ggsave("../figures/strauss_quant1.pdf", plot = s2.plot, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)

t_max = 0.25
l = which(L_df$t ==t_max)
s3.plot = ggplot(data = s.quantiles2[1:l,]) +
  geom_line(aes(y = lower, x = L_df$t[1:l]), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t[1:l]), col = 'black')+
  geom_line(data = L_cells_df[1:l,], aes(y = L, x = t), col = 'red' )+
  ggtitle("Cells L-function with Strauss event RF 0.01 and 0.99 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

s3.plot

ggsave("../figures/strauss_quant2.pdf", plot = s3.plot, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 8/3, units = "in",
       dpi = 300, limitsize = TRUE)


