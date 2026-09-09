## return posterior predictive plots for the genes given the parameter fits
## return a list of figs along each axis


posterior_predictive_plots <- function(gene, y, par, r_x, r_t, r_tx,
                                       x_observed, x_test, P_obs,
                                       dsqx, dsqt, dsqtx, B, N, M,
                                       slides_per_group, noise_diag, 
                                       kernel = 'RBF', 
                                       fit_donor = F,
                                       cov_donor,
                                       meta, 
                                       slide_df, 
                                       FILTER = FALSE # whether to filter to contrast
                                       ){
                                       
  ## function arguements are
  # par : GP parameters for the gene / gene
  # obs / test coordinates
  # SE kernel pairwise distance (dsqx, dsqt)
  # B : coregionalization matrix
  # noise_diag : slide-specific noise  
  # kernel : RBF / Matern --- what cov kernel to fit to data
  
  ## extract GP parameters
  bw = par[1]
  sigma_f = par[2]
  sigma_d = ifelse(fit_donor, par[3], 0)
  sigma_n = par[4:length(par)]
  n_slides <- sum(slides_per_group)
  Sigma <- diag(sigma_n[noise_diag] ^ 2)
  
  ## compute covariance kernels
  if(kernel == 'RBF'){
    Kxx <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2))
    Ktt <- sigma_f^2 * exp(-0.5 * dsqt / (bw^2))
    Ktx <- sigma_f^2 * exp(-0.5* dsqtx / (bw^2))
  }else if(kernel == 'Matern5'){
    ## compute explicitly for each case. only two computations for each kernel
    Kxx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
    Ktt <- sigma_f^2 * (1 + (sqrt(5)*r_t / bw) + 5*dsqt/(3*bw^2)) * exp(-(sqrt(5)*r_t / bw))
    Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_tx / bw) + 5*dsqtx/(3*bw^2)) * exp(-(sqrt(5)*r_tx / bw))
  } else if(kernel == 'Matern3'){
    Kxx <- sigma_f^2 * (1 + sqrt(3)*r_x / bw) * exp(-sqrt(3)*r_x / bw)
    Ktt <- sigma_f^2 * (1 + sqrt(3)*r_t / bw) * exp(-sqrt(3)*r_t / bw)
    Ktx <- sigma_f^2 * (1 + sqrt(3)*r_tx / bw) * exp(-sqrt(3)*r_tx / bw)
  }else{
    stop(sprintf("Unsupported kernel type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", kernel))
  }
  
  
  ## covariance for observations
  Kf <- kronecker(B, Kxx) + kronecker(Sigma, diag(N)) 
  Kf <- P_obs %*% tcrossprod(Kf, P_obs) + sigma_d^2*cov_donor
  U_f <- chol(Kf)
  Kf_inv <- chol2inv(U_f)
  
  ##--------------------------------------------------
  ## posterior mean for observations (UW^T | Y)
  ##--------------------------------------------------
  ## if we want donor-specific shifted curves for each mouse
  if(fit_donor){
    donor_levels <- unique(slide_df$donor_id)
    test_donors <- rep(slide_df$donor_id, each=M)
    Z_test <- model.matrix(~0 + factor(test_donors, levels = donor_levels))
    Z_obs <- model.matrix(~0 + factor(meta$donor_id, levels = donor_levels))
    
    cov_donor_star <- tcrossprod(Z_test, Z_obs)
    Ktx_obs <- tcrossprod(kronecker(B, Ktx) , P_obs) + sigma_d^2*cov_donor_star
  }else{
    Ktx_obs <- tcrossprod(kronecker(B, Ktx) , P_obs)
  }
  ##--------------------------------------------------
  ## posterior predictives
  y_pred <- Ktx_obs %*% Kf_inv %*% y

  ##--------------------------------------------------
  ##plot the data
  ##--------------------------------------------------
  gp_df <- data.frame(
    y_gp    = as.vector(y_pred),
    coord_x = rep(x_test, times = n_slides),
    slide   = rep(slide_df$batch, each = M),
    donor = rep(slide_df$donor_id, each=M),
    age     = rep(slide_df$age, each = M)
  )

  ## filter to only time points of interest
  ## subsetting after i have the posterior fits
  if(FILTER){
    gp_df <- filter(gp_df, age %in% c('24wk', '90wk'))
    meta <- filter(meta, age %in% c('24wk', '90wk'))
  }
  
  
  # 2. Compute the mean trajectory across slides within each age group
  gp_mean_df <- gp_df %>%
    group_by(age, coord_x) %>%
    summarize(y_gp = mean(y_gp), .groups = "drop")

  # 3. Plot the clean mean curve
  p <- ggplot(data = meta, aes(x = coord_x, y = !!sym(gene), fill = donor_id)) +
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
      aes(x = coord_x, y = y_gp, group = age),
      color = 'red', linewidth = 0.75, alpha = .85,
      inherit.aes=F
    ) + 
    ## posterior fits for each donor
    # geom_line(
    #   data = gp_df,
    #   aes(x = coord_x, y = y_gp, colour = donor, group = slide),
    #   linewidth = 0.75, alpha = .85, linetype = 1,
    #   inherit.aes = FALSE
    # ) +
    facet_wrap(~ age) +
    ylim(quantile(meta %>% pull(gene), c(0.005, 0.995), na.rm = TRUE)) +
    labs(x = 'Isodepth', y = gene, fill = "Donor") +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 16, face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(size = 16)
    )
  
  return(p)
}
