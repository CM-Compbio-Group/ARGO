plot_posterior <- function(y_pred, argo_env){
  
  with(argo_env, {
    
  # if(FIT_DONOR){
  #   donor_levels <- levels(metadata$donor_id)
  #   test_donors <- rep(metadata %>% distinct(sample_id, donor_id) %>% pull(donor_id), each=M)
  #   Z_test <- model.matrix(~0 + factor(test_donors, levels = donor_levels))
  #   Z_obs <- model.matrix(~0 + metadata$donor_id)
  #   cov_donor_star <- tcrossprod(Z_test, Z_obs)
  #   Ktx_obs <- tcrossprod(kronecker(B, Ktx) , projection) + sigma_d^2*cov_donor_star
  # }
  ##--------------------------------------------------
  ##plot the data
  ##--------------------------------------------------
  slide_df <- metadata %>% distinct(sample_id, donor_id, condition)
  grid_df <- expand.grid(
    coord_x   = test_coords,
    sample_id = slide_df$sample_id
  )
  grid_df <- merge(grid_df, slide_df, by = "sample_id", sort = FALSE)
  
  gp_df <- data.frame(
    y_gp    = as.vector(y_pred),
    coord_x = grid_df$coord_x,
    slide   = grid_df$sample_id,
    donor = grid_df$donor_id,
    condition = grid_df$condition
  )
  
  ## filter to only time points of interest
  ## subsetting after i have the posterior fits
  plot_meta <- bind_cols(metadata, feature = feature_vec)
  
  gp_df <- filter(gp_df, condition %in% c('24wk', '90wk'))
  plot_meta <- filter(plot_meta, condition %in% c('24wk', '90wk'))
  
  
  # 2. Compute the mean trajectory across slides within each condition group
  gp_mean_df <- gp_df %>%
    group_by(condition, coord_x) %>%
    summarize(y_gp = mean(y_gp), .groups = "drop")
  
  
  # 3. Plot the clean mean curve
  p <- ggplot(data = plot_meta, 
              aes(x = coord_x, y = feature, color = donor_id)) +
    geom_jitter(
      width = 0.015,
      height = 0,
      shape = 16,
      alpha = 0.35,
      size = .75
    ) 
  
  ## posterior fits for each donor
  if(FIT_DONOR){
    p <- p + geom_line(
      data = gp_df,
      aes(x = coord_x, y = y_gp, color = donor, group = slide),
      linewidth = 0.75, alpha = .85, linetype = 1,
      inherit.aes = FALSE
    )
  }
  
  p <- p + # Single mean GP curve per facet
    geom_line(
      data = gp_mean_df,
      aes(x = coord_x, y = y_gp, group = condition),
      color = 'red', linewidth = 0.5, alpha = .75,
      inherit.aes=F
    ) + 
    facet_wrap(~ condition) +
    ylim(quantile(feature_vec, c(0.005, 0.995), na.rm = TRUE)) +
    labs(x = 'Isodepth', y = .feature, color = "Donor") +
    scale_x_continuous(breaks = c(0, 0.5, 1.0), 
                       labels = c("0.0", "0.5", "1.0")) + 
    theme_bw(base_size = 8, 
             base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), # remove grid lines
      
      legend.box.spacing = unit(1, "mm"),           # Drastically reduces space between plot border and legend
      legend.margin = margin(l = -3, r = 0, t = 0, b = 0, unit = "pt"), # Pulls the legend closer to the left
      legend.key.size = unit(3, "mm"), # Optional: compacts color keys vertically too
      
      # 2. UNIFORM 8 PT PLAIN FONTS EVERYWHERE
      text = element_text(size = 8, face = "plain"),
      plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
      axis.title = element_text(size = 8, face = "plain"),
      axis.text = element_text(size = 6, face = "plain"),
      
      # Facet labels (e.g., "24wk", "90wk", and the gene names on the side)
      strip.text = element_text(size = 8, face = "plain"), 
      # strip.background = element_blank(),           # Optional: Cleans up facet header boxes
      
      # Legend text and title formatting
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8, face = "plain")
      
    )
  
  return(p)
  })
}