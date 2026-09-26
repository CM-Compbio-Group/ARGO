library(ggplot2)

# Parameters from your setup
amplitude <- 1
bandwidth <- 0.1 # Try changing this to 0.025
center <- 0.5
N_bins <- 20

# Simulate pseudotime coordinates (t) from 0 to 1
t_vals <- seq(0, 1, length.out = 1000)
lfc_vals <- amplitude * exp(-((t_vals - center)^2) / (2 * bandwidth^2))

# Define bin edges
bin_edges <- seq(0, 1, length.out = N_bins + 1)

# Plot
df_plot <- data.frame(t = t_vals, lfc = lfc_vals)
p <- ggplot(df_plot, aes(x = t, y = lfc)) +
  geom_line(color = "blue", size = 1.2) +
  geom_vline(xintercept = bin_edges, linetype = "dashed", color = "gray50", alpha = 0.6) +
  labs(
    title = glue::glue("Perturbation Profile (Bandwidth = {bandwidth}, Bins = {N_bins})"),
    subtitle = "Dotted lines represent pseudotime bin boundaries",
    x = "Pseudotime (t)",
    y = "Log-Fold Change (LFC)"
  ) +
  theme_minimal()

plot(p)