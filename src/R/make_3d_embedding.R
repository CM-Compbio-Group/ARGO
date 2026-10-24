fx <- 'cOXPHOS'

stat <- funGP_fits[[fx]]$statistic
band <- funGP_fits[[fx]]$band
df <- tibble(bin = levels(bin_meta$bin), 
             sig_domain =  ifelse((stat > band) | (stat < -band), 1, 0)) %>%
  mutate(sigDomainAlpha = ifelse(sig_domain == 1, 0.9, .1),
         sigDomainSize = ifelse(sig_domain == 1, 5, 5))

bin_WT <- left_join(bin_WT, df,by = 'bin')
bin_Mutant <- left_join(bin_Mutant, df, by = 'bin')

all_colors <- c(bin_Mutant[[fx]], bin_WT[[fx]])
c_min <- min(all_colors, na.rm = TRUE)
c_max <- max(all_colors, na.rm = TRUE)

color_palette <- scales::col_numeric(palette = "viridis", domain = c(c_min, c_max))
bin_WT$hex_color <- color_palette(bin_WT[[fx]])
bin_Mutant$hex_color <- color_palette(bin_Mutant[[fx]])

x_range <- range(c(bin_WT$CTmvp, bin_Mutant$CTmvp), na.rm = TRUE)
y_range <- range(c(bin_WT$CTpan, bin_Mutant$CTpan), na.rm = TRUE)
z_range <- range(c(bin_WT$LE, bin_Mutant$LE), na.rm = TRUE)

p_wt <- plotly::plot_ly(type = "scatter3d", mode = "markers") %>%
  # Significant points
  plotly::add_markers(
    data = dplyr::filter(bin_WT, sig_domain == 1),
    x = ~CTmvp, y = ~CTpan, z = ~LE, 
    showlegend = FALSE,
    scene = "scene2", # Target the second 3D viewport
    marker = list(
      size = 5, 
      opacity = 0.5, 
      color=~hex_color
      # color = as.formula(paste0("~", fx)), # Keep mapping formula
      # coloraxis = "coloraxis"               # Bind to single layout scale
    )
  ) %>%
  # Non-significant points
  plotly::add_markers(
    data = dplyr::filter(bin_WT, sig_domain == 0),
    x = ~CTmvp, y = ~CTpan, z = ~LE, 
    showlegend = FALSE,
    scene = "scene2",
    marker = list(
      size = 2.5, 
      opacity = 0.15, 
      color=~hex_color
      # color = as.formula(paste0("~", fx)),
      # coloraxis = "coloraxis"
    )
  )


p_mutant <- plotly::plot_ly(type = "scatter3d", mode = "markers") %>%
  plotly::add_markers(
    data = dplyr::filter(bin_Mutant, sig_domain == 1),
    x = ~CTmvp, y = ~CTpan, z = ~LE, 
    showlegend = FALSE,
    scene = "scene", # Target the first 3D viewport
    marker = list(
      size = 5, 
      opacity = 0.5, 
      color=~hex_color
      # color = as.formula(paste0("~", fx)),
      # coloraxis = "coloraxis"
    )
  ) %>%
  # Non-significant points
  plotly::add_markers(
    data = dplyr::filter(bin_Mutant, sig_domain == 0),
    x = ~CTmvp, y = ~CTpan, z = ~LE, 
    showlegend = FALSE,
    scene = "scene",
    marker = list(
      size = 2.5, 
      opacity = 0.15, 
      color=~hex_color
      # color = as.formula(paste0("~", fx)),
      # coloraxis = "coloraxis"
    )
  ) 


# -------------------------------------------------------------------------
# 2. Build Standardized Structural Settings for Scenes
# -------------------------------------------------------------------------
shared_scene_settings <- list(
  aspectmode = "cube",  # Grid shapes are identical boxes
  xaxis = list(title = "CTmvp", range = x_range),
  yaxis = list(title = "CTpan", range = y_range),
  zaxis = list(title = "LE", range = z_range)
)

final_plot <- plotly::subplot(p_mutant, p_wt, nrows = 1, margin = 0.05) %>%
  plotly::layout(
    title = paste("Abundance analysis -", fx),
    scene = c(list(domain = list(x = c(0, 0.45))), shared_scene_settings),
    scene2 = c(list(domain = list(x = c(0.55, 1.0))), shared_scene_settings)
  )

final_plot$x$layout <- final_plot$x$layout[grep('NA', names(final_plot$x$layout), invert = TRUE)]


final_plot
