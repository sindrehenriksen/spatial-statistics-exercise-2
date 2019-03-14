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

# Plot prior realizations
prior_sim_plots = vector("list", n_sims)
margin = theme(plot.margin = unit(c(0.4, 0, -0.4, 0), "cm"))
labs = labs(x="", y="")
d = data
for(i in 1:n_sims){
  d["pines"] = prior_sims[,i]
  prior_sim_plots[[i]] = ggplot(d, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill + margin + labs
}
prior_sim_grid_plot = arrangeGrob(grobs=prior_sim_plots, ncol=2)
ggsave("../figures/p2_prior_sims.pdf", plot=prior_sim_grid_plot,
       width=5, height=8, units="in", dpi=300)

# Plot simulated observations####
prior_sim_obs_plots = vector("list", n_sims)
for(i in 1:n_sims){
  d["pines"] = rbinom(n, prior_sims[,i], data$alpha)
  prior_sim_obs_plots[[i]] = ggplot(d, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill
}

# Generate 10 realizations from posterior with estimated intensity
posterior_intensity = node_area*lambda_hat*(1 - data$alpha)
posterior_sims = replicate(n_sims, rpois(n, posterior_intensity))

# Plot posterior realizations
posterior_sim_plots = vector("list", n_sims)
for(i in 1:n_sims){
  d["pines"] = posterior_sims[,i] + data$pines
  posterior_sim_plots[[i]] = ggplot(d, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill + margin + labs
}
posterior_sim_grid_plot = arrangeGrob(grobs=posterior_sim_plots, ncol=2)
ggsave("../figures/p2_posterior_sims.pdf", plot=posterior_sim_grid_plot,
       width=5, height=8, units="in", dpi=300)

# Plot simulated observations####
posterior_sim_obs_plots = vector("list", n_sims)
for(i in 1:n_sims){
  d["pines"] = rbinom(n, posterior_sims[,i], data$alpha)
  posterior_sim_obs_plots[[i]] = ggplot(d, aes(x, y)) +
    geom_raster(aes(fill=pines)) + fill
}

# Plot posterior expected number of pine trees
d["pines"] = posterior_intensity + data$pines
posterior_mean_plot = ggplot(data=d, aes(x, y)) +
  geom_raster(aes(fill=pines)) + fill
