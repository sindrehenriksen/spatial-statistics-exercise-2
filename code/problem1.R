#Spatial statistics
#Problem 1
###############################################################
###############################################################

#Load libraries
library(spatial)
library(ggplot2)
library(gridExtra)

#Read data
#NB these are a bit different from the ones on the web page. 
#y has different sign in redwood
#pines does not look to be the same data

#cells = ppinit('cells.dat')
#redwood = ppinit('redwood.dat')
#pines = ppinit('pines.dat')

#Data from web page 
#Not all on the same form as above..
cells = as.list(read.table('../data_files/cells.dat', col.names = c('x', 'y')))
redwood = as.list(read.table('../data_files/redwood.dat', col.names = c('x', 'y')))
#pines was on a different form..
pines = ppinit('../data_files/pines.dat')
###############################################################
###############################################################
#a)

#Display cells data
cells_df = data.frame(x = cells$x, y = cells$y)

ggplot(data = cells_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Cells data") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figures/prob1_cells_points.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
      dpi = 300, limitsize = TRUE)

#Display redwood data
redwood_df = data.frame(x = redwood$x, y = redwood$y)

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figures/prob1_redwood_points.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

#Display pines data
pines_df = data.frame(x = pines$x, y = pines$y)

ggplot(data = pines_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Pines data") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figures/prob1_pines_points.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 4, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

###############################################################
###############################################################
#b)
#Kfn function: Computes L = sqrt(K/pi). Cannot see how this corresponds quite to J(t)?
L_plot = list()

L_cells = Kfn(cells, fs = 1)
J_cells = L_cells$y/L_cells$x
L_cells_df = data.frame(t = L_cells$x, L = L_cells$y, J = J_cells)

L_plot[[1]] = ggplot(data = L_cells_df, aes(x=t, y=L)) + 
  geom_path() + 
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle("L-interaction function for cells") +
  theme(plot.title = element_text(hjust = 0.5))


L_redwood = Kfn(redwood, fs = 1) 
J_redwood = L_redwood$y/L_redwood$x
L_redwood_df = data.frame(t = L_redwood$x, L = L_redwood$y, J = J_redwood)

L_plot[[2]] = ggplot(data = L_redwood_df, aes(x=t, y=L)) + 
  geom_path() + 
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle("L-interaction function for redwood trees") +
  theme(plot.title = element_text(hjust = 0.5))


L_pines = Kfn(pines, fs = 1)
J_pines = L_pines$y/L_pines$x
L_pines_df = data.frame(t = L_pines$x, L = L_pines$y, J = J_pines)

L_plot[[3]] = ggplot(data = L_pines_df, aes(x=t, y=L)) + 
  geom_path() + 
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle("L-interaction function for pine trees") +
  theme(plot.title = element_text(hjust = 0.5))

L_emp = grid.arrange(grobs = L_plot, ncol = 1)

ggsave("../figures/prob1_L_empirical.pdf", plot = L_emp, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4*2, units = "in",
       dpi = 300, limitsize = TRUE)

###############################################################
###############################################################
#c)

