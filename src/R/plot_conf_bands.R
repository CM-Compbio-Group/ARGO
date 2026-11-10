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
    # p <- ggplot(plot_data, aes(x = x)) +
    #   geom_line(aes(y = y), linewidth = 0.5) +
    #   geom_line(aes(y = band_upper), color = "red", linetype = "dashed", linewidth = 0.5) +
    #   geom_line(aes(y = band_lower), color = "red", linetype = "dashed", linewidth = 0.5) +
    #   coord_cartesian(ylim = y_limits) +
    #   labs(
    #     x = "Isodepth",
    #     y = bquote(Delta * " expression"),
    #     title = .feature
    #     # title = glue("{conf_level_alpha*100}% Confidence band for {.feature}")
    #   ) +
    #   theme_bw(base_size = 8, 
    #                 base_family = "Arial") + 
    #   theme( 
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     
    #     # Forces the plot title to be exactly 8pt, plain (not bold), and centered
    #     plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    #     
    #     # Forces axis labels ("Isodepth", "LFC") to be exactly 8pt and plain
    #     axis.title = element_text(size = 8, face = "plain"),
    #     
    #     # Forces axis tick numbers (0.00, 0.25, etc.) to be exactly 8pt and plain
    #     axis.text = element_text(size = 8, face = "plain")
    #   )
    
    p <- ggplot(plot_data, aes(x = x)) +
      geom_ribbon(aes(ymin = band_lower, ymax = band_upper, fill = "Confidence Band"), alpha = 0.2) +
      geom_line(aes(y = y, color = "Test statistic"), linewidth = 0.5, linetype = 2) +
      scale_fill_manual(name = NULL, values = c("Confidence Band" = "red3")) +
      scale_color_manual(name = NULL, values = c("Test statistic" = "black")) +
      coord_cartesian(ylim = y_limits) +
      labs(
        x = "Isodepth",
        y = bquote(Delta * " expression"),
        title = .feature
      ) +
      theme_bw(base_size = 8, base_family = "Arial") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        # Uses a negative top margin to pull the legend up toward the plot
        legend.box.margin = margin(t = -8, r = 0, b = 0, l = 0, unit = "pt"),
        plot.margin = margin(t = 0.5, r = 2, b = 0.5, l = 2, unit = "pt")
      )
    
    
    return(p)
  })
}
    