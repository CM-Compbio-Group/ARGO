plot_CI_bands <- function(gene, x, band, stat, alpha){
  # 1. Combine your data into a single dataframe
  plot_data <- data.frame(
    x = x,
    y = stat,
    band_upper = band,
    band_lower = -band
  )
  
  # 2. Dynamically calculate y-limits
  y_limits <- range(c(plot_data$y, plot_data$band_upper, plot_data$band_lower)) 
  
  # 3. Create the ggplot2 chart
  p <- ggplot(plot_data, aes(x = x)) +
    geom_line(aes(y = y), linewidth = 0.8) +
    geom_line(aes(y = band_upper), color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_line(aes(y = band_lower), color = "red", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(ylim = y_limits) +
    labs(
      x = "Isodepth value",
      y = "Mean Difference",
      title = glue("{alpha}% Confidence band for {gene}")
    ) +
    theme_minimal()
  
  return(p)
}