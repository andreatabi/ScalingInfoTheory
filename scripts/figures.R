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
load("~/sims_block.RData")


## Figure S2 =============================================================================================

DKL_curve_over_alpha_k <- function(
    N0, gamma, phi,
    alpha_grid,
    k,
    nsim = 1e5,
    nbins = 100,
    disc = "equalfreq",
    smooth_span = 0.25
) {
  df_r <- purrr::map_dfr(alpha_grid, function(a) {
    sim <- adult_B_samples_raw(Nmax = N0, alpha = a, nsim = nsim, gamma = 1, phi = phi)
  #  rk  <- coarsegrain_block_sum(sim$r, k)
    Rk  <- coarsegrain_block_sum(sim$r, k)
    rk  <- rescale_clt(Rk,k)
    tibble::tibble(alpha = a, r = rk)
  })
  
  r_disc <- infotheo::discretize(df_r$r, disc = disc, nbins = nbins)[, 1]
  alpha_disc <- as.factor(df_r$alpha)
  
  out <- local_MI_alpha_discrete(alpha_disc, r_disc, logbase = 2) %>% dplyr::mutate(Nmax = N0, gamma = gamma, phi = phi,
                                                                       k = k, nsim =nsim, nbins = nbins)
  
  if (!is.null(smooth_span) && length(unique(out$alpha)) > 4) {
    fit <- stats::loess(D_KL ~ alpha, data = out, span = smooth_span)
    out$D_KL_smooth <- stats::predict(fit, newdata = data.frame(alpha = out$alpha))
  } else {
    out$D_KL_smooth <- out$D_KL
  }
  
  out
}

alpha_grid <- seq(0.9, 2, length.out=100)
phi <- 10^(-5.3) # beta ~ 0.75
k_grid <- c(1,2,4,8,16)

dk_curves <- purrr::map_dfr(k_grid, ~ local_KL_cg(
  N0 = 1e14, gamma = 1, phi = phi,
  alpha_grid = alpha_grid,
  k = .x,
  nsim = 2e5,
  nbins = 100
))
head(dk_curves)

ggplot(dk_curves, aes(alpha, D_KL_smooth, group = factor(k), linetype = factor(k))) +
  geom_line() +
  theme_classic() +
  labs(linetype = "k (block size)", y = "KL-to-mixture (smoothed)")

alpha_mins <- map_dfr(k_grid, ~find_alpha_min_cg(N0=1e14, gamma=1, 
                                phi=10^(-5.3), 
                                k= .x, 
                                alpha_grid=alpha_grid, 
                                nsim=2e5, nbins=100))
alpha_mins

ks <- c(1,2,4,8,16)
sds <- sapply(ks, function(k){
  sim <- adult_B_samples_raw(Nmax=1e14, alpha=1.79, nsim=2e5, gamma=1, phi=10^(-5.3))
  Rk <- coarsegrain_block_sum(sim$r, k)
  sd(Rk)
})
data.frame(k=ks, sd=sds, sd_over_sqrtk=sds/sqrt(ks))

ggplot(dk_curves,aes(x = alpha,y = D_KL_smooth, group = factor(k),linetype = factor(k)) )+
  geom_line() +
  geom_point(data=alpha_mins, aes(alpha_min, KL_min),
             size = 2,
             color = "black",
             fill  = "white",
              shape = 21,
             stroke = 0.7) +
  scale_color_brewer(palette = "Spectral",name= "Body size (log10)") +
  labs(x = expression(alpha), 
       y = expression( D[KL](alpha, r ) ),
       linetype = "Block size (k)"  ) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  theme(legend.position   = c(0.8,0.7) ,panel.grid.minor  = element_blank() )+
  theme(aspect.ratio = 1)

ggsave("~/FigS2.png", width = 4, height = 4)

## Figure 2A ===========================================================================
N_grid <- c(1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16)
alpha_grid <- seq(0.9, 2, length.out=100)

beta_mins <- map_dfr(N_grid, ~ find_alpha_beta_min2(N0=.x, gamma=1, 
                                                    phi=10^(-5.3), 
                                                    alpha_grid=alpha_grid, 
                                                    nsim=1e5, nbins=100))
beta_mins
mean(beta_mins$beta_min)

df_curves <- purrr::map_dfr(N_grid, ~ local_MI_alpha_r_for_species_params(
  N0 = .x, gamma = 1, 
  phi = 10^(-5.3),
  alpha_grid = alpha_grid,
  nsim = 1e5,
  nbins = 100
))
head(df_curves)

ggplot(df_curves,aes(x = alpha,y = D_KL_smooth, color = factor(log10(Nmax)),group = factor(log10(Nmax))) )+
  geom_line() +
  geom_point(data=beta_mins, aes(alpha_min, D_KL_min),
             size = 2,
             color = "black",
             fill  = "white",
             shape = 21,
             stroke = 0.7) +
  scale_color_brewer(palette = "Spectral",name= "Body size (log10)") +
  labs(x = expression(alpha), 
       y = expression( D[KL](alpha, r ) )) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  theme(legend.position   = "none" ,panel.grid.minor  = element_blank() )+
  theme(aspect.ratio = 1)

ggsave("~/Fig2a.png", width = 4, height = 4)

## Figure 2B ===========================================================================

df_smooth <- df_beta_min_surface %>%
  group_by(Nmax) %>%
  group_modify(~{
    dd <- .
    dd$logphi <- log10(dd$phi)
    fit_a <- loess(alpha_min ~ logphi, data = dd, span = 0.4)
    fit_b <- loess(beta_min ~ logphi, data = dd, span = 0.4)
    dd$alpha_smooth <- predict(fit_a, newdata = dd)
    dd$beta_smooth <- predict(fit_b, newdata = dd)
    dd
  }) %>%
  ungroup()

df_deriv <- df_smooth %>%
  arrange(Nmax, logphi) %>%
  group_by(Nmax) %>%
  mutate(
    d_alpha_dlogphi = c(NA, diff(alpha_smooth) / diff(logphi) )
  ) %>%
  ungroup()
df_deriv

delta <- 0.01

phi_bands_species <- df_deriv %>%
  group_by(Nmax) %>%
  filter(!is.na(d_alpha_dlogphi)) %>%
  filter(abs(d_alpha_dlogphi) < delta) %>%  # nearly flat
  summarise(
    phi_min   = min(phi),
    phi_max   = max(phi),
    alpha_star_max = max(alpha_smooth),
    alpha_star = mean(alpha_smooth),
    beta_star  = mean(beta_smooth, na.rm = TRUE),
    beta_sd  = sd(beta_smooth, na.rm = TRUE),
    n = n(),
    .groups   = "drop"
  )
phi_bands_species

phi_star_per_species <- df_smooth %>%
    group_by(Nmax) %>%
    filter(!is.na(alpha_min)) %>%
    slice_max(order_by = alpha_smooth, with_ties = TRUE) %>%
    ungroup() %>%
    select(Nmax, phi_star = phi, alpha_star = alpha_smooth,
           D_KL_min, beta_star = beta_smooth)



ggplot(df_smooth, aes(x = phi, y = alpha_smooth, color = factor(Nmax))) +
  geom_line(size = 0.75) +
  geom_point( data = phi_star_per_species, aes(x = phi_star, y = alpha_star),
              inherit.aes = FALSE, size = 2, color = "black", fill="white", shape = 21, stroke = 0.7 ) +
  scale_x_log10() +
  labs(x = expression(phi),
       y = expression(alpha[min]),
       color = "Body size" ) +
  scale_color_brewer(palette = "Spectral") +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position  = "none",
    aspect.ratio     = 1
  )

ggsave("~/fig2b.png", width = 4, height = 4)


# Figure 2C =================================================================
df_univ <- df_smooth %>%
  group_by(Nmax) %>%
  mutate(
    alpha_star = (alpha_smooth - min(alpha_smooth, na.rm = TRUE)) /
      (max(alpha_smooth, na.rm = TRUE) - min(alpha_smooth, na.rm = TRUE)),
    
    beta_star  = (beta_smooth  - min(beta_smooth,  na.rm = TRUE)) /
      (max(beta_smooth,  na.rm = TRUE) - min(beta_smooth,  na.rm = TRUE))
  ) %>%
  ungroup()


ggplot(df_univ,aes(x = alpha_star, y = beta_star, 
                   group = factor(log10(Nmax)), color = factor(log10(Nmax)))) +
  geom_path( size=0.75) +
  scale_color_brewer(palette="Spectral")+
  labs(x = expression(alpha[min]^"'"), y = expression(beta^"'" * (alpha[min])),
       color = expression("Species size " * (log[10]))) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  theme(aspect.ratio = 1)

ggsave("~/fig2c.png", width = 4, height = 4)

# Figure 2D =================================================================

load("~/simsFI.RData")

df_reduced <- df_reduced %>%
  filter(is.finite(FI_min), FI_min > 0,
         is.finite(MI_min))

df_reduced_beta <- df_reduced %>%
  mutate(
    beta_min = beta_species(
      Nmax  = Nmax,
      alpha = alpha_min,
      gamma = 1,
      phi   = phi
    )
  ) %>%
  filter(is.finite(MI_min),
         is.finite(FI_min),
         is.finite(beta_min),
         FI_min > 0)


df_reduced %>% filter(alpha_min == 1.7)

library(wesanderson)

ggplot(df_reduced_beta, aes(x = FI_min, y = MI_min, color = beta_min) ) +
  geom_point(pch=21, size=2, alpha = 0.9) +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))+
  labs(
    x     = expression(I[r](alpha[min]) ),
    y     = expression(D[KL](alpha[min], r) ),
    shape = expression("Body size" * (log[10]) ),
    color = expression(beta(alpha[min])) ) +
  theme(legend.position   = "right",
        panel.grid.minor  = element_blank(),
        plot.title        = element_text(face = "bold"))+
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        legend.position = c(0.85, 0.4),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  theme(aspect.ratio = 1)


ggsave("~/Fig2d.png", width = 4, height = 4)


### Figure S1 =================================================================
phi_grid <- 10^seq(-6, -2, length.out = 100) 
N_grid     <- c(1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16)
alpha_grid <- seq(0.9, 2.0, length.out = 100)


beta_SOGM <- function(alpha, N_grid) {
  out <- lapply(N_grid, function(Nmax) {
    beta <- beta_species(Nmax, alpha, gamma = 1, phi=1)
  })
  out
}

dfbeta <- expand_grid(
  Nmax  = N_grid,
  alpha = alpha_grid,
  phi   = phi_grid
) %>%
  mutate(
    beta = mapply(
      function(N, a, p) beta_species(
        Nmax  = N,
        alpha = a,
        gamma = 1,
        phi   = p
      ),
      Nmax, alpha, phi
    )
  )

dfbeta


ggplot(dfbeta, aes(alpha, beta, color=(phi), group=factor(phi))) +
  geom_line() +
  facet_wrap(~log10(Nmax))+
  scale_color_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 1), linetype="dashed")+
  labs(x = expression(alpha), y = expression(beta[S]), color = expression(phi)) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position = c(0.5,0.15))


ggsave("~/FigS1.png", width = 6, height = 7)




