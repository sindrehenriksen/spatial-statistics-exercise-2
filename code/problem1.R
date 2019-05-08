#Spatial statistics
#Problem 1
###############################################################
###############################################################

#Load libraries
library(spatial)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

#Read data

#cells = ppinit('cells.dat')
#redwood = ppinit('redwood.dat')
#pines = ppinit('pines.dat')

#NB these are a bit different from the ones on the web page. 
#y has different sign in redwood
#pines does not look to be the same data

#Data from web page 
#Not all on the same form as above..
cells = as.list(read.table('../data/cells.dat', col.names = c('x', 'y')))
redwood = as.list(read.table('../data/redwood.dat', col.names = c('x', 'y')))
#pines was on a different form..
pines = ppinit('../data/pines.dat')
###############################################################
###############################################################
#a)

#Display cells data
cells_df = data.frame(x = cells$x, y = cells$y)

ggplot(data = cells_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Cells data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#ggsave("../figures/prob1_cells_points.pdf", plot = last_plot(), device = NULL, path = NULL,
#       scale = 1, width = 4, height = 4, units = "in",
#      dpi = 300, limitsize = TRUE)

#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("../figures/prob1_redwood_points.pdf", plot = last_plot(), device = NULL, path = NULL,
#       scale = 1, width = 4, height = 4, units = "in",
#      dpi = 300, limitsize = TRUE)

#Display pines data
pines_df = data.frame(x = pines$x, y = pines$y)

ggplot(data = pines_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Pines data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("../figures/prob1_pines_points.pdf", plot = last_plot(), device = NULL, path = NULL,
#       scale = 1, width = 4, height = 4, units = "in",
#       dpi = 300, limitsize = TRUE)

###############################################################
###############################################################
#b)
#Kfn function: Computes L = sqrt(K/pi). Cannot see how this corresponds quite to J(t)?
L_plot = list()

L_cells = Kfn(cells, fs = 1)
J_cells = L_cells$y/L_cells$x
L_cells_df = data.frame(t = L_cells$x, L = L_cells$y, J = J_cells, 
                        Data = rep('cells', length(J_cells)))

L_plot[[1]] = ggplot(data = L_cells_df, aes(x=t, y=L)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_path(col = 'red') + 
  ggtitle("L-interaction function for cells") +
  theme(plot.title = element_text(hjust = 0.5))


L_redwood = Kfn(redwood, fs = 1) 
J_redwood = L_redwood$y/L_redwood$x
L_redwood_df = data.frame(t = L_redwood$x, L = L_redwood$y, J = J_redwood,
                          Data = rep('redwood', length(J_cells)))

L_plot[[2]] = ggplot(data = L_redwood_df, aes(x=t, y=L)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_path(col = 'red') + 
  ggtitle("L-interaction function for redwood trees") +
  theme(plot.title = element_text(hjust = 0.5))


L_pines = Kfn(pines, fs = 1)
J_pines = L_pines$y/L_pines$x
L_pines_df = data.frame(t = L_pines$x, L = L_pines$y, J = J_pines,
                        Data = rep('pines', length(J_cells)))

L_plot[[3]] = ggplot(data = L_pines_df, aes(x=t, y=L)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_path(col = 'red') + 
  ggtitle("L-interaction function for pine trees") +
  theme(plot.title = element_text(hjust = 0.5))

#Empirical and theoretical L
L_emp_theor = grid.arrange(grobs = L_plot, ncol = 1)

#ggsave("../figures/prob1_L_emp_theor.pdf", plot = L_emp_theor, device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 4*2, units = "in",
#       dpi = 300, limitsize = TRUE)

#All empirical in one display
L_all = rbind(L_cells_df, L_redwood_df, L_pines_df)

ggplot(data = L_all, aes(x=t, y=L, col = Data)) + 
  geom_path() + 
  ggtitle("L-interaction functions") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("../figures/prob1_L_empirical.pdf", plot = last_plot(), device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 4, units = "in",
#       dpi = 300, limitsize = TRUE)
###############################################################
###############################################################
#c)
set.seed(4250)

#Condition on number of observed points
Num_cells = length(cells$x)
Num_redw = length(redwood$x)
Num_pines = length(pines$x)
Num_samps = 100

#Relevant positions in both directions are [0,1]
#Generate realizations
samps_cells_x = matrix(runif(Num_samps*Num_cells), ncol = Num_samps)
samps_cells_y = matrix(runif(Num_samps*Num_cells), ncol = Num_samps)

samps_redwood_x = matrix(runif(Num_samps*Num_redw), ncol = Num_samps)
samps_redwood_y = matrix(runif(Num_samps*Num_redw), ncol = Num_samps)

samps_pines_x = matrix(runif(Num_samps*Num_pines), ncol = Num_samps)
samps_pines_y = matrix(runif(Num_samps*Num_pines), ncol = Num_samps)


#Preallocate matrix to save L-function
L_samps_cells_t = matrix(0, nrow = 70, ncol = Num_samps)
L_samps_cells_L = matrix(0, nrow = 70, ncol = Num_samps)

L_samps_redwood_t = matrix(0, nrow = 70, ncol = Num_samps)
L_samps_redwood_L = matrix(0, nrow = 70, ncol = Num_samps)

L_samps_pines_t = matrix(0, nrow = 70, ncol = Num_samps)
L_samps_pines_L = matrix(0, nrow = 70, ncol = Num_samps)

#Compute L for each sample
for (i in 1:Num_samps){
  temp_L_cells = Kfn(list(x = samps_cells_x[,i], 
                          y = samps_cells_y[,i]), fs = 1)
  L_samps_cells_t[,i] = temp_L_cells$x
  L_samps_cells_L[,i] = temp_L_cells$y
  
  temp_L_redwood = Kfn(list(x = samps_redwood_x[,i], 
                            y = samps_redwood_y[,i]), fs = 1)
  L_samps_redwood_t[,i] = temp_L_redwood$x
  L_samps_redwood_L[,i] = temp_L_redwood$y
  
  temp_L_pines = Kfn(list(x = samps_pines_x[,i], 
                          y = samps_pines_y[,i]), fs = 1)
  L_samps_pines_t[,i] = temp_L_pines$x
  L_samps_pines_L[,i] = temp_L_pines$y
}

#############################################
#Empirical 0.9-intervals
cells.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
redwood.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
pines.quantiles = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  cells.quantiles[i,]= quantile(L_samps_cells_L[i,], c(0.05,0.95))
  redwood.quantiles[i,]= quantile(L_samps_cells_L[i,], c(0.05,0.95))
  pines.quantiles[i,]= quantile(L_samps_cells_L[i,], c(0.05,0.95))
}

#Empirical 0.98-intervals
cells.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
redwood.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
pines.quantiles2 = data.frame(lower = rep(0,70), upper =rep(0,70))
for (i in 1:70){
  cells.quantiles2[i,]= quantile(L_samps_cells_L[i,], c(0.01,0.99))
  redwood.quantiles2[i,]= quantile(L_samps_cells_L[i,], c(0.01,0.99))
  pines.quantiles2[i,]= quantile(L_samps_cells_L[i,], c(0.01,0.99))
}

#############################################
#Displaying:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Cells:
L_df = as.data.frame(L_samps_cells_L)
L_df$t = L_samps_cells_t[,1]
#Changing format to be able to plot all realisations
long_L = melt(L_df, id = 't')

cells1 = ggplot(long_L,
             aes(x=t, y=value, colour=variable)) +
  scale_colour_manual(values=rep(cbPalette, length = Num_samps))+
  geom_line()+
  geom_line(data = cells.quantiles, aes(y = lower, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  geom_line(data = cells.quantiles, aes(y = upper, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  #geom_line(data = L_cells_df, aes(y = L), col = 'black' )+
  ggtitle("Generated L-functions: Cells") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

#melt(data.frame(L_cells_df$t, L_cells_df$L, cells.quantiles), id = 't')
cells2 = ggplot(data = cells.quantiles) +
  geom_line(aes(y = lower, x = L_df$t), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t), col = 'black')+
  geom_line(data = L_cells_df, aes(y = L, x = t), col = 'red' )+
  ggtitle("Cells L-function with Poisson RF 0.05 and 0.95 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))
###########################
#Redwood:
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

###########################
#Pines:
L_df = as.data.frame(L_samps_pines_L)
L_df$t = L_samps_pines_t[,1]
#Changing format to be able to plot all realisations
long_L = melt(L_df, id = 't')

pines1 = ggplot(long_L,
       aes(x=t, y=value, colour=variable)) +
  scale_colour_manual(values=rep(cbPalette, length = Num_samps))+
  geom_line()+
  geom_line(data = pines.quantiles, aes(y = lower, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  geom_line(data = pines.quantiles, aes(y = upper, x = L_df$t), 
            inherit.aes = FALSE, col = 'black')+
  #geom_line(data = L_pines_df, aes(y = L), col = 'black' )+
  ggtitle("Generated L-functions: Pines") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

pines2 = ggplot(data = pines.quantiles) +
  geom_line(aes(y = lower, x = L_df$t), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t), col = 'black')+
  geom_line(data = L_pines_df, aes(y = L, x = t), col = 'red' )+
  ggtitle("Pines L-function with Poisson RF 0.05 and 0.95 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

plot.samples = grid.arrange(grobs = list(cells1, redwood1, pines1), ncol = 1)

#ggsave("../figures/prob1_samples.pdf", plot = plot.samples, device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 4*2, units = "in",
#       dpi = 300, limitsize = TRUE)

plot.quantiles = grid.arrange(grobs = list(cells2, redwood2, pines2), ncol = 1)

ggsave("../figures/prob1_quantiles.pdf", plot = plot.quantiles, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4*2, units = "in",
      dpi = 300, limitsize = TRUE)


#98% intervals: 
t_max = 0.25
l = which(L_df$t ==t_max) 
cells3 = ggplot(data = cells.quantiles2[1:l,]) +
  geom_line(aes(y = lower, x = L_df$t[1:l]), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t[1:l]), col = 'black')+
  geom_line(data = L_cells_df[1:l,], aes(y = L, x = t), col = 'red' )+
  ggtitle("Cells L-function with Poisson RF 0.01 and 0.99 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

redwood3 = ggplot(data = redwood.quantiles2[1:l,]) +
  geom_line(aes(y = lower, x = L_df$t[1:l]), col = 'black')+
  geom_line(aes(y = upper, x = L_df$t[1:l]), col = 'black')+
  geom_line(data = L_redwood_df[1:l,], aes(y = L, x = t), col = 'red' )+
  ggtitle("Redwood L-function with Poisson RF 0.01 and 0.99 quantiles") +
  xlab("t")+
  ylab("L") +
  theme(plot.title = element_text(hjust = 0.5))

plot.quantiles2 = grid.arrange(grobs = list(cells3, redwood3), ncol = 1)

#ggsave("../figures/prob1_quantiles2.pdf", plot = plot.quantiles2, device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 4*2, units = "in",
#      dpi = 300, limitsize = TRUE)