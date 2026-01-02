rm(list=ls())
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

#### SWEEP  ==========================================================
library(progressr)
handlers(global = TRUE)   
handlers("txtprogressbar")  

phi_grid <- 10^seq(-6, -2, length.out = 100) 
N_grid     <- c(1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16)
alpha_grid <- seq(0.9, 2.0, length.out = 100)

param_grid <- expand.grid(
  Nmax  = N_grid,
  gamma = 1,
  phi   = phi_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% as_tibble()

df_beta_min_surface <- with_progress({
  
  p <- progressor(steps = nrow(param_grid)) 
  
  pmap_dfr(param_grid,function(Nmax, gamma, phi) {
    p()
    
    find_alpha_beta_min2(
      N0         = Nmax,
      gamma      = gamma,
      phi        = phi,
      alpha_grid = alpha_grid,
      nsim       = 1e5,   
      nbins      = 150
    )
  })
})

save(df_beta_min_surface, file = "~/sims.RData")





