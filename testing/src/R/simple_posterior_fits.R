## return posterior predictive plots for the genes given the parameter fits
simple_posterior_fit <- function(gene) {
  
  ## pick out variable and parameters
  y <- bin_meta %>% pull(gene)
  par <- params[gene, ]
  ## extract GP parameters
  bw = par[1]
  sigma_f = par[2]
  sigma_n = par[3:length(par)]
  
  K <- length(levels(bin_meta$condition))
  n_slides <- n_samples * K
  Sigma <- diag(sigma_n[noise_diag] ^ 2)
  
  ## compute covariance KERNELs
  if(KERNEL == 'RBF'){
    Kxx <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2))
    Ktt <- sigma_f^2 * exp(-0.5 * dsqt / (bw^2))
    Ktx <- sigma_f^2 * exp(-0.5* dsqtx / (bw^2))
  }else if(KERNEL == 'Matern5'){
    ## compute explicitly for each case. only two computations for each KERNEL
    Kxx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
    Ktt <- sigma_f^2 * (1 + (sqrt(5)*r_t / bw) + 5*dsqt/(3*bw^2)) * exp(-(sqrt(5)*r_t / bw))
    Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_tx / bw) + 5*dsqtx/(3*bw^2)) * exp(-(sqrt(5)*r_tx / bw))
  } else if(KERNEL == 'Matern3'){
    Kxx <- sigma_f^2 * (1 + sqrt(3)*r_x / bw) * exp(-sqrt(3)*r_x / bw)
    Ktt <- sigma_f^2 * (1 + sqrt(3)*r_t / bw) * exp(-sqrt(3)*r_t / bw)
    Ktx <- sigma_f^2 * (1 + sqrt(3)*r_tx / bw) * exp(-sqrt(3)*r_tx / bw)
  }else{
    stop(sprintf("Unsupported KERNEL type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", KERNEL))
  }
  
  
  ## covariance for observations
  Kf <- kronecker(B, Kxx) + kronecker(Sigma, diag(N)) 
  U_f <- chol(Kf)
  Kf_inv <- chol2inv(U_f)
  ##--------------------------------------------------
  ## posterior mean for observations (UW^T | Y)
  ##--------------------------------------------------
  Ktx_obs <- kronecker(B, Ktx)
  ##--------------------------------------------------
  ## posterior predictives
  y_pred <- Ktx_obs %*% Kf_inv %*% y
  ##--------------------------------------------------
  ##plot the data
  ##--------------------------------------------------
  gp_df <- data.frame(
    y_gp    = as.vector(y_pred),
    coord_x = rep(x_test, times = n_slides),
    condition     = rep(levels(bin_meta$condition), each = n_samples * M)
  )
  
  # 2. Compute the mean trajectory across slides within each age group
  gp_mean_df <- gp_df %>%
    group_by(condition, coord_x) %>%
    summarize(y_gp = mean(y_gp), .groups = "drop")
  
  # 3. Plot the clean mean curve
  p <- ggplot(data = bin_meta, aes(x = coord_x, y = !!sym(gene), fill = sample)) +
    geom_jitter(
      width = 0.015,
      height = 0,
      shape = 21,
      color = "transparent",
      alpha = 0.5,
      size = 1.5
    ) +
    # Single mean GP curve per facet
    geom_line(
      data = gp_mean_df,
      aes(x = coord_x, y = y_gp, group = condition),
      color = 'red', linewidth = 0.75, alpha = .85,
      inherit.aes=F
    ) + 
    # ggtitle(gene) + 
    facet_wrap(~ condition) +
    ylim(quantile(bin_meta %>% pull(gene), c(0.005, 0.995), na.rm = TRUE)) +
    labs(x = 'Isodepth', y = gene, fill = "sample") +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 16, face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(size = 16)
    )
  
  return(p)
}