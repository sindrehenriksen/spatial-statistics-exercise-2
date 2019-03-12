setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(tidyverse)
# library(spatial)
# library(MASS)

# data = ppinit("pines.dat")
data = as_tibble(read.table("../data/obsprob.txt", header=TRUE))
data["pines"] = read.table("../data/obspines.txt", header=TRUE)[3]
pine_coords = tibble(x=rep(data$x, data$pines), y=rep(data$y, data$pines))

p = ggplot(data, aes(x, y))
fill = scale_fill_distiller(palette='Spectral')
pines_plot = p + geom_raster(aes(fill=pines)) + fill
alpha_plot = p + geom_raster(aes(fill=alpha)) + fill
ggsave("../figures/p1_pines_plot.pdf", plot=pines_plot,
       width=4, height=4, units="in", dpi=300)
ggsave("../figures/p1_alpha_plot.pdf", plot=alpha_plot,
       width=4, height=4, units="in", dpi=300)
