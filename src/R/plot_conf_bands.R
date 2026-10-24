plot_CI_bands <- function(feature, res, argo_env){
  argo_env$.feature <- feature
  argo_env$.stat <- res$statistic
  argo_env$.band <- res$band
  with(argo_env, {
    
    plot_data <- data.frame(
      x = test_coords,
      y = .stat,
      band_upper = .band,
      band_lower = -.band
    )
    
    # 2. Dynamically calculate y-limits
    y_limits <- range(c(plot_data$y, plot_data$band_upper, plot_data$band_lower)) 
    
    # 3. Create the ggplot2 chart
    p <- ggplot(plot_data, aes(x = x)) +
      geom_line(aes(y = y), linewidth = 0.5) +
      geom_line(aes(y = band_upper), color = "red", linetype = "dashed", linewidth = 0.5) +
      geom_line(aes(y = band_lower), color = "red", linetype = "dashed", linewidth = 0.5) +
      coord_cartesian(ylim = y_limits) +
      labs(
        x = "Isodepth",
        y = "LFC",
        title = .feature
        # title = glue("{conf_level_alpha*100}% Confidence band for {.feature}")
      ) +
      theme_bw(base_size = 8, 
                    base_family = "Arial") + 
      theme( 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        # Forces the plot title to be exactly 8pt, plain (not bold), and centered
        plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
        
        # Forces axis labels ("Isodepth", "LFC") to be exactly 8pt and plain
        axis.title = element_text(size = 8, face = "plain"),
        
        # Forces axis tick numbers (0.00, 0.25, etc.) to be exactly 8pt and plain
        axis.text = element_text(size = 8, face = "plain")
      )
    
    
    return(p)
  })
}
    