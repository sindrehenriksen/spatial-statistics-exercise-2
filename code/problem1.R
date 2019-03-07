#Spatial statistics
#Problem 1
###############################################################
###############################################################

#Load libraries
library(spatial)
library(ggplot2)

#Read data
cells = ppinit('cells.dat')
redwood = ppinit('redwood.dat')
pines = ppinit('pines.dat')


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

#Er nok ikke helt riktig
L_cells = Kfn(cells, fs = 1)
#L_redwood = Kfn(redwood, )
L_pines = Kfn(pines, fs = 10)


###############################################################
###############################################################
#c)

