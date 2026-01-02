rm(list=ls())
library(hrbrthemes)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(extraDistr)
library(moments)
library(VGAM)
library(infotheo)
library(tibble)
source("~/functions.R")
load("~/sims.RData")


### FISHER INFO ========================================================
N_grid     <- c(1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16)
alpha_grid <- seq(0.9, 2, length.out = 200)

df_fi <- map_dfr(N_grid, function(N0) {
    fisher_r_alpha_discrete(
      N0  = N0,
      phi = 10^(-5.3),
      alpha_grid = alpha_grid,
      gamma   = 1,
      nsim    = 2e5,
      nbins = 100
    )
})

beta_mins <- map_dfr(N_grid, ~ find_alpha_beta_min2(N0=.x, gamma=1, 
                                                    phi=10^(-5.3), 
                                                    alpha_grid=alpha_grid, 
                                                    nsim=2e5, nbins=100))
beta_mins
min(df_fi$FI_r, na.rm = T)

ggplot(df_fi, aes(x = alpha,y = FI_r, color = factor(log10(Nmax)),
                 group = factor(log10(Nmax))) )+
  geom_line() +
  scale_color_brewer(palette = "Spectral",name=expression("Species size " * (log[10]))) +
  labs(x = expression(alpha), 
       y = expression( I[r](alpha) )) +
  geom_vline(data= beta_mins, aes(xintercept=alpha_min, color = factor(log10(Nmax))), linetype="dashed")+
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  theme(panel.grid.minor  = element_blank() )+
  theme(aspect.ratio = 1)

ggsave("~/Dropbox/Scaling_InfoTheory/figures/FigS3.png", width = 5, height = 4)


N_grid     <- c(1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16)
alpha_grid <- seq(0.9, 2, length.out = 100)
phi_grid   <- 10^seq(-6, -2, length.out = 100)

df_reduced <- map_dfr(N_grid, function(N0) {
  map_dfr(phi_grid, function(phi_val) {
    compute_MI_FI_min_for_phi(
      N0         = N0,
      phi_val    = phi_val,
      alpha_grid = alpha_grid,
      gamma      = 1,
      nsim_MI    = 1e5,
      nsim_FI    = 1e5,
      nbins_MI   = 100,
      nbins_FI   = 100
    )
  })
})

save(df_reduced, file = "~/Dropbox/Scaling_InfoTheory/data/simsFI.RData")

