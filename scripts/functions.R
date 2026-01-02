

B_adult_alpha <- function(Nmax, Bc, Ec, gamma, phi, eta, dt = 0.1, kneg = 0.3, alpha = 1) {
  s   <- phi * (Nmax)^(alpha/2) * sqrt(dt / 2)
  add <- s * ( Ec * ( 1/(2*eta) + (1 - eta + kneg)/2 ) + gamma )
  (Nmax * Bc + add) / 3600
}

m_adult <- function(Nmax, Bc, Ec, phi, m_c, eta, dt = 0.1) {
  s <- phi * sqrt(Nmax) * sqrt(dt / 2)
  m_c * ( Nmax + (Ec / (Bc * eta)) * (s / 2) )
}

beta_species <- function(Nmax, alpha, Bc = 3e-11 * 3600, Ec = 3.7e-7, 
                         gamma = 1, phi = 1, eta = 0.7, dt = 0.1, kneg = 0.3) {
  A  <- Bc
  C0 <- phi * sqrt(dt / 2) * ( Ec * ( 1/(2*eta) + (1 - eta + kneg)/2 ) + gamma )
  E0 <- phi * sqrt(dt / 2) * ( (Ec / Bc) * (1/(2*eta)) )
  n <- Nmax
  
  num1 <- A * n + (alpha / 2) * C0 * n^(alpha / 2)
  den1 <- A * n + C0 * n^(alpha / 2)
  num2 <- n + E0 * n^(alpha / 2)
  den2 <- n + (alpha / 2) * E0 * n^(alpha / 2)
  beta <- (num1 / den1) * (num2 / den2)
  beta
}

adult_B_samples <- function(Nmax, alpha = 1, nsim = 2e5, Bc = 3e-11 * 3600, Ec = 3.7e-7, 
                            gamma = 1, phi= 1, 
                            eta=0.7, dt = 0.1, kneg = 0.3) {
  B0 <- Nmax * Bc
  s  <- phi * (Nmax)^(alpha/2) * sqrt(dt / 2)
  
  x  <- VGAM::rlaplace(nsim, 0, s)
  P  <- pmax(x, 0)  
  N  <- pmax(-x, 0)  
  
  c_pos <- Ec/eta + (1 - eta) * Ec + gamma
  c_neg <- kneg * Ec + gamma
  
  B  <- B0 + c_pos * P + c_neg * N
  B  <- pmax(B, 1e-20) 

  tau <- 1
  r   <- log(B[(1 + tau):length(B)] / B[1:(length(B) - tau)])
  r_lap <- mean(abs(r - median(r)))
  r_var <- var(r)
  kurt  <- kurtosis(r)
  
  list(mean_logB = mean(log(B)),
       mean_B = mean(B),
       mean_B_det = mean(B0),
       var_mean  = var(log(B)) / nsim,
       r_mean = mean(r),
       r_var = r_var, 
       kurt = kurt ) 
}


adult_B_samples_raw <- function(Nmax, alpha, nsim = 2e5,
                                gamma, phi, 
                                Bc = 3e-11 * 3600, Ec = 3.7e-7, 
                                eta = 0.7, dt = 0.1, kneg = 0.3) {
  B0 <- Nmax * Bc
  s  <- phi * (Nmax)^(alpha/2) * sqrt(dt / 2)
  
  x  <- VGAM::rlaplace(nsim, 0, s)
  P  <- pmax(x, 0)
  N  <- pmax(-x, 0)
  
  c_pos <- Ec/eta + (1 - eta) * Ec + gamma
  c_neg <- kneg * Ec + gamma
  
  B  <- B0 + c_pos * P + c_neg * N
  B  <- pmax(B, 1e-20) 
  
  tau <- 1L
  if (length(B) > tau) {
    r <- log(B[(1 + tau):length(B)] / B[1:(length(B) - tau)])
  } else {
    r <- numeric(0)
  }
  
  list(B = B, r = r)
}


beta_species <- function(Nmax, alpha, gamma, phi, Bc = 3e-11 * 3600, Ec = 3.7e-7, 
                         eta=0.7, dt = 0.1, kneg = 0.3) {
  A  <- Bc
  C0 <- phi * sqrt(dt / 2) * ( Ec * ( 1/(2*eta) + (1 - eta + kneg)/2 ) + gamma )
  E0 <- phi * sqrt(dt / 2) * ( (Ec / Bc) * (1/(2*eta)) )
  
  n <- Nmax
  num <- (A + C0 / (2 * n^(alpha/2)  )) * (n + E0 * n^(alpha/2) )
  den <- (A * n + C0 *  n^(alpha/2) ) * (1 + E0 / (2 * n^(alpha/2) ))
  num / den
}


local_MI_alpha_discrete <- function(alpha_disc, x_disc, logbase = 2) {
  tab   <- table(alpha_disc, x_disc)
  joint <- tab / sum(tab)
  
  p_alpha <- rowSums(joint)
  p_x     <- colSums(joint)
  
  D_KL_vec <- apply(joint, 1, function(row_joint) {
    if (sum(row_joint) == 0) return(NA_real_)
    p_x_given_a <- row_joint / sum(row_joint)
    idx <- p_x_given_a > 0 & p_x > 0
    sum(p_x_given_a[idx] * (
      log(p_x_given_a[idx] / p_x[idx]) / log(logbase)
    ))
  })
  
  MI_global <- sum(p_alpha * D_KL_vec, na.rm = TRUE)
  
  tibble(
    alpha     = as.numeric(rownames(joint)),
    p_alpha   = as.numeric(p_alpha),
    D_KL      = as.numeric(D_KL_vec),
    contrib   = p_alpha * D_KL_vec,
    MI_global = MI_global
  )
}

local_MI_alpha_r_for_species_params <- function(
    N0, alpha_grid,
    gamma, phi,
    nsim = 2e4,
    nbins = 100,
    smooth_span = 0.25  
) {
  df_r <- map_dfr(alpha_grid, function(a) {
    sim <- adult_B_samples_raw(
      Nmax  = N0,
      alpha = a,
      nsim  = nsim,
      gamma = gamma,
      phi   = phi
    )
    tibble(alpha = a, r = sim$r)
  })
  
  r_disc     <- discretize(df_r$r, disc = "equalfreq", nbins = nbins)[, 1]
  alpha_disc <- as.factor(df_r$alpha)
  
  out <- local_MI_alpha_discrete(alpha_disc, r_disc, logbase = 2) %>%
    mutate(Nmax = N0, gamma = gamma, phi = phi)
  
  if (!is.null(smooth_span)) {
    out <- out %>%
      arrange(alpha)
    
    fit <- loess(D_KL ~ alpha, data = out, span = smooth_span)
    #fit <- gam(D_KL ~ s(alpha, k = 10), data = out, method = "REML")
    
    out$D_KL_smooth <- as.numeric(predict(fit, newdata = data.frame(alpha = out$alpha)))
  }
  
  out
}

find_alpha_beta_min2 <- function(
    N0, gamma=1, phi,
    alpha_grid,
    nsim = 2e4,
    nbins = 100,
    smooth_span = 0.25   
) {
  df_local <- local_MI_alpha_r_for_species_params(
    N0         = N0,
    alpha_grid = alpha_grid,
    gamma      = gamma,
    phi        = phi,
    nsim       = nsim,
    nbins      = nbins
  )
  
  df_local <- df_local %>%
    arrange(alpha)
  
  if (length(unique(df_local$alpha)) > 4) {
    fit <- loess(D_KL ~ alpha, data = df_local, span = smooth_span)
    df_local$D_KL_smooth <- predict(fit, newdata = df_local)
    
    na_idx <- is.na(df_local$D_KL_smooth)
    df_local$D_KL_smooth[na_idx] <- df_local$D_KL[na_idx]
    
    row_min <- df_local %>%
      slice_min(order_by = D_KL_smooth, n = 1, with_ties = FALSE)
    
    MI_min   <- row_min$D_KL_smooth
  } else {
    row_min <- df_local %>%
      slice_min(order_by = D_KL, n = 1, with_ties = FALSE)
    
    MI_min   <- row_min$D_KL
  }
  
  alpha_min <- row_min$alpha
  
  beta_min  <- beta_species(
    Nmax  = N0,
    alpha = alpha_min,
    gamma = gamma,
    phi   = phi
  )
  
  tibble(
    Nmax      = N0,
    gamma     = gamma,
    phi       = phi,
    alpha_min = alpha_min,
    D_KL_min    = MI_min,
    beta_min  = beta_min,
    MI_global = unique(df_local$MI_global),
    nbins = nbins,
    nsim = nsim
  )
}

fisher_r_alpha_discrete <- function(
    N0, alpha_grid,
    gamma, phi,
    nsim  = 2e4,
    nbins = 100
) {
  df_r <- map_dfr(alpha_grid, function(a) {
    sim <- adult_B_samples_raw(
      Nmax  = N0,
      alpha = a,
      nsim  = nsim,
      gamma = gamma,
      phi   = phi
    )
    tibble(alpha = a, r = sim$r)
  })
  
  r_disc     <- discretize(df_r$r, disc = "equalfreq", nbins = nbins)[, 1]
  alpha_disc <- factor(df_r$alpha, levels = alpha_grid)
  
  tab   <- table(alpha_disc, r_disc)
  joint <- tab / sum(tab)             # p(alpha, r_bin)
  
  # Conditional p(r_bin | alpha) as matrix
  p_alpha <- rowSums(joint)
  cond    <- joint / p_alpha          # rows sum to 1
  cond[cond <= 0] <- 1e-12            # avoid zeros
  
  # 3. Fisher via curvature of log p(r|alpha)
  n_alpha <- length(alpha_grid)
  delta   <- mean(diff(alpha_grid))
  
  FI_vec <- rep(NA_real_, n_alpha)
  for (i in 2:(n_alpha - 1)) {
    p_i   <- cond[i,  ]                 # p(x | alpha_i)
    lp_p  <- log(cond[i + 1, ])        # log p(x | alpha_{i+1})
    lp_m  <- log(cond[i - 1, ])        # log p(x | alpha_{i-1})
    dlogp <- (lp_p - lp_m) / (2 * delta)
    
    FI_vec[i] <- sum(p_i * dlogp^2)
  }
  
  tibble(
    Nmax  = N0,
    gamma = gamma,
    phi   = phi,
    alpha = alpha_grid,
    FI_r  = FI_vec
  )
}


compute_MI_FI_min_for_phi <- function(
    N0, phi_val,
    alpha_grid,
    gamma,
    nsim_MI  = 1e5,
    nsim_FI  = 1e5,
    nbins_MI = 100,
    nbins_FI = 100,
    smooth_span = 0.25
) {
  df_MI <- local_MI_alpha_r_for_species_params(
    N0         = N0,
    alpha_grid = alpha_grid,
    gamma      = gamma,
    phi        = phi_val,
    nsim       = nsim_MI,
    nbins      = nbins_MI
  )
  
  df_MI <- df_MI %>%
    arrange(alpha) %>%
    mutate(
      D_KL_smooth = {
        fit <- loess(D_KL ~ alpha, span = smooth_span)
        predict(fit, newdata = data.frame(alpha = alpha))
      }
    )
  
  row_min <- df_MI %>%
    slice_min(order_by = D_KL_smooth, n = 1, with_ties = FALSE)
  
  alpha_min <- row_min$alpha
  MI_min    <- row_min$D_KL_smooth
  
  df_FI <- fisher_r_alpha_discrete(
    N0         = N0,
    alpha_grid = alpha_grid,
    gamma      = gamma,
    phi        = phi_val,
    nsim       = nsim_FI,
    nbins      = nbins_FI
  )
  
  FI_min <- approx(
    x = df_FI$alpha,
    y = df_FI$FI_r,
    xout = alpha_min
  )$y
  
  tibble(
    Nmax      = N0,
    phi       = phi_val,
    gamma     = gamma,
    alpha_min = alpha_min,
    MI_min    = MI_min,
    FI_min    = FI_min
  )
}


coarsegrain_block_mean <- function(r, k) {
  n <- length(r)
  m <- floor(n / k)
  if (m < 2) return(numeric(0))
  r2 <- r[1:(m*k)]
  rowMeans(matrix(r2, nrow = k, ncol = m))
}

coarsegrain_block_sum <- function(r, k) {
  n <- length(r); m <- floor(n / k)
  if (m < 2) return(numeric(0))
  r2 <- r[1:(m*k)]
  colSums(matrix(r2, nrow = k, ncol = m))
}

rescale_clt <- function(Rk, k) Rk / sqrt(k)

rescale_z <- function(x) as.numeric(scale(x))

local_KL_cg <- function(
    N0, alpha_grid, gamma, phi,
    nsim = 2e4, nbins = 100,
    k = 1,
    smooth_span = 0.25
) {
  df_rk <- purrr::map_dfr(alpha_grid, function(a) {
    sim <- adult_B_samples_raw(Nmax = N0, alpha = a, nsim = nsim, gamma = gamma, phi = phi)
    #rk  <- coarsegrain_block_sum(sim$r, k)
    Rk <- coarsegrain_block_sum(sim$r, k)
    #rk <- rescale_z(Rk)        
    rk <- rescale_clt(Rk,k)        
    tibble::tibble(alpha = a, r = rk)
  })
  
  r_disc <- infotheo::discretize(df_rk$r, disc = "equalfreq", nbins = nbins)[, 1]
  alpha_disc <- factor(df_rk$alpha, levels = alpha_grid)
  
  out <- local_MI_alpha_discrete(alpha_disc, r_disc, logbase = 2) %>%
    dplyr::mutate(Nmax = N0, gamma = gamma, phi = phi, k = k) %>%
    dplyr::arrange(alpha)
  
  if (!is.null(smooth_span) && length(unique(out$alpha)) > 4) {
    fit <- stats::loess(D_KL ~ alpha, data = out, span = smooth_span)
    out$D_KL_smooth <- stats::predict(fit, newdata = data.frame(alpha = out$alpha))
    out$D_KL_smooth[is.na(out$D_KL_smooth)] <- out$D_KL[is.na(out$D_KL_smooth)]
  } else {
    out$D_KL_smooth <- out$D_KL
  }
  
  out
}


find_alpha_min_cg <- function(
    N0, gamma, phi, alpha_grid,
    nsim = 2e4, nbins = 100,
    k = 1, smooth_span = 0.25
) {
  df <- local_KL_cg(
    N0 = N0, alpha_grid = alpha_grid, gamma = gamma, phi = phi,
    nsim = nsim, nbins = nbins, k = k, smooth_span = smooth_span
  )
  
  row_min <- df %>% dplyr::slice_min(order_by = D_KL_smooth, n = 1, with_ties = FALSE)
  
  beta_min  <- beta_species(
    Nmax  = N0,
    alpha = row_min$alpha,
    gamma = gamma,
    phi   = phi
  )
  
  tibble::tibble(
    Nmax = N0, gamma = gamma, phi = phi, k = k,
    alpha_min = row_min$alpha,
    KL_min = row_min$D_KL_smooth,
    beta = beta_min,
    nsim = nsim, nbins = nbins
  )
}
