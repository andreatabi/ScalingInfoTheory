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
k_grid <-  c(1,2,4,8,16)

param_grid <- expand.grid(
  Nmax  = N_grid,
  gamma = 1,
  phi   = phi_grid,
  k = k_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% as_tibble()

df <- with_progress({
  
  p <- progressor(steps = nrow(param_grid))
  
  pmap_dfr(param_grid, function(Nmax, gamma, phi, k) {
    p()
    
    find_alpha_min_cg(
      N0         = Nmax,
      gamma      = gamma,
      phi        = phi,
      k          = k,
      alpha_grid = alpha_grid,
      nsim       = 2e5,
      nbins      = 100
    )
  })
})

head(df)

save(df, file = "~/sims_block.RData")


#### test ============================================================================================
nbins_for_k <- function(k) max(20, floor(100 / sqrt(k)))

test1 <- find_alpha_min_cg(N0=1e12, gamma=1, phi=1e-5, k=1, alpha_grid=alpha_grid, nsim=1e5, nbins=100)
test2 <- find_alpha_min_cg(N0=1e12, gamma=1, phi=1e-5, k=8, alpha_grid=alpha_grid, nsim=1e5, nbins=100)
test3 <- find_alpha_min_cg(N0=1e12, gamma=1, phi=1e-5, k=16, alpha_grid=alpha_grid, nsim=1e5, nbins=100)

rbind(test1, test2, test3)


summarize_valley <- function(df_kl, window = 0.08) {
  df_kl <- dplyr::arrange(df_kl, alpha)
  y <- df_kl$D_KL_smooth %||% df_kl$D_KL
  
  i0 <- which.min(y)
  alpha0 <- df_kl$alpha[i0]
  
  Delta <- max(y) - min(y)
  
  local <- dplyr::filter(df_kl, abs(alpha - alpha0) <= window)
  if (nrow(local) >= 5) {
    fit <- stats::lm(D_KL_smooth ~ poly(alpha, 2, raw = TRUE), data = local)
    b2 <- stats::coef(fit)[3]  # coefficient of alpha^2
    curvature <- 2 * b2
  } else {
    curvature <- NA_real_
  }
  
  tibble::tibble(alpha_min = alpha0, Delta = Delta, curvature = curvature)
}

summarize_valley(df)

