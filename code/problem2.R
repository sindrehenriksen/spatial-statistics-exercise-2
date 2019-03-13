setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(tidyverse)
# library(spatial)
# library(MASS)

# Load data
data = as_tibble(read.table("../data/obsprob.txt", header=TRUE))
data["pines"] = read.table("../data/obspines.txt", header=TRUE)[3]
# pine_coords = tibble(x=rep(data$x, data$pines), y=rep(data$y, data$pines))
tot_area = 300^2
node_area = 100
n = nrow(data)

# Plot number of observed pines and observation probabilities
p = ggplot(data, aes(x, y))
fill = scale_fill_distiller(palette='Spectral')
pines_plot = p + geom_raster(aes(fill=pines)) + fill
alpha_plot = p + geom_raster(aes(fill=alpha)) + fill
ggsave("../figures/p2_pines.pdf", plot=pines_plot,
       width=5, height=4, units="in", dpi=300)
ggsave("../figures/p2_alpha.pdf", plot=alpha_plot,
       width=5, height=4, units="in", dpi=300)

# Calculate estimate of intensity
lambda_hat = sum(data$pines / data$alpha) / tot_area

# Generate 10 realizations from prior with estimated intensity
n_sims = 10
prior_sims = replicate(n_sims, rpois(n, node_area*lambda_hat))

# Plot realizations
prior_sim_plots = vector("list", n_sims)
margin = theme(plot.margin = unit(c(0.4, 0, -0.4, 0), "cm"))
labs = labs(x="", y="")
for(i in 1:n_sims){
  d = tibble(x=data$x, y=data$y, pines=prior_sims[,i])
  prior_sim_plots[[i]] = ggplot(d, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill + margin + labs
}
prior_sim_plots[[1]] = prior_sim_plots[[1]] + top_margin
prior_sim_plots[[2]] = prior_sim_plots[[2]] + top_margin
prior_sim_grid_plot = arrangeGrob(grobs=prior_sim_plots, ncol=2)
ggsave("../figures/p2_prior_sims.pdf", plot=prior_sim_grid_plot,
       width=5, height=8, units="in", dpi=300)

# Plot simulated observations
sim_obs_plots = vector("list", n_sims)
for(i in 1:n_sims){
  o = tibble(x=data$x, y=data$y, pines=rbinom(n, prior_sims[,i], data$alpha))
  prior_sim_obs_plots[[i]] = ggplot(o, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill
}
