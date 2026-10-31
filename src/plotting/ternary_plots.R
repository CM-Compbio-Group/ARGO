library(ggpubr)
library(cowplot)
library(gridExtra)
library(extrafont)

make_tern_plot <- function(df, condition, state){
  
  df$inROI <- ifelse(df$bin %in% roi_index[[state]], 'T', 'F')
  p <- ggtern(data = df, aes(x = simplex_CTpan, 
                                    y = simplex_CTmvp, 
                                    z = simplex_LE, 
                                    color = !!sym(state), 
                                    size = inROI, 
                                    alpha = inROI)) + 
    labs("T" = "CTmvp", 
         "L" = "CTpan",
         "R" = "LE", 
         x=NULL,y=NULL,z=NULL,
         title = condition) + 
    geom_point() +
    theme_nogrid_major() + 
    scale_color_gradientn(
      colors = c("lightgray", "#fde0dd", "#fa9fb5", "#c51b8a", "#50115a"),
      limits = c(-1, 2),
      name = state,
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "right",
        barheight = unit(1, "in"),
        barwidth = unit(0.08, "in"),
        frame.colour = "black",
        order = 1
      )
    ) + 
    scale_alpha_manual(values = c('F' = .25, 'T' = 0.75)) + 
    scale_size_manual(values = c('F' = 0.5, 'T' = 3)) + 
    theme(tern.panel.background = element_blank(), 
          text = element_text(family = "Arial", size = 8)
          )
  return(p)
}


ivy_tern <- c('CTpan', 'CTmvp', 'LE')

# 1. Compute bin centroids from the UNCLIPPED raw module scores in brain@meta.data
#    so negative scores retain their continuous gradient instead of a flat ReLU floor
raw_bin_coords <- brain@meta.data %>%
  mutate(bin = factor(km$cluster)) %>% 
  group_by(bin) %>%
  summarise(across(all_of(ivy_tern), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# 2. Convert to [0, 1] relative scores per axis (Min-Max or Empirical CDF quantiles)
#    Raising to power `gamma` (e.g., 1.5 - 2) controls how strongly points push toward edges/vertices
gamma <- 1.5
tern_mat <- raw_bin_coords %>%
  select(all_of(ivy_tern)) %>%
  as.matrix()

# Min-max scale each of the 3 axes to [0, 1] across bins so the minimum touches the triangle edge (0)
tern_scaled <- apply(tern_mat, 2, function(x) {
  ((x - min(x)) / (max(x) - min(x)) + 0.02)^gamma
})

# Normalize rows to sum to 1 on the 2-simplex (triangle)
simplex_coords <- sweep(tern_scaled, 1, rowSums(tern_scaled), '/')
colnames(simplex_coords) <- paste('simplex', colnames(simplex_coords), sep = '_')

simplex_df <- bind_cols(bin = raw_bin_coords$bin, as.data.frame(simplex_coords))

# 3. Join onto bin_meta (replacing any previous simplex_ columns)
bin_meta <- bin_meta %>%
  select(-any_of(c("simplex_CT", "simplex_CTmvp", "simplex_CTpan", "simplex_LE"))) %>%
  left_join(simplex_df, by = "bin")


roi_index <- lapply(funGP_fits, function(kxx) which(abs(kxx$statistic) > kxx$band))
sapply(roi_index, length)

stats <- funGP_fits$cMESInf$statistic
band <- funGP_fits$cMESInf$band
idx <- which((stats > band) | (stats < -band) )

plot(simplex_df$simplex_CTpan, stats, ylim = range(c(stats, band, -band)), col = 'blue', pch=20, cex=1.5)
points(simplex_df$simplex_CTpan, band, col = 'red', pch = 20)
points(simplex_df$simplex_CTpan, -band, col = 'red', pch = 20)



# 4. Plot on ggtern
bin_meta$inROI <- ifelse(bin_meta$bin %in% roi_index$cMESInf, 'T', 'F')


state <- 'cMESVasc'
p1 <- make_tern_plot(filter(bin_meta, condition == 'IDHmutant'), 'IDHmutant', state)
p2 <- make_tern_plot(filter(bin_meta, condition == 'GBM'), 'IDHwt', state)

p1_clean <- p1 + theme(legend.position = "none", plot.margin = margin(15, 10, 10, 25))
p2_clean <- p2 + theme(legend.position = "none", plot.margin = margin(15, 10, 10, 10))

shared_legend <- get_legend(p1)
g1 <- ggplotGrob(p1_clean)
g2 <- ggplotGrob(p2_clean)

# 4. Use patchwork's wrap_elements() to arrange the frozen grobs and legend side-by-side
p_final <- wrap_elements(full = g1) + 
  wrap_elements(full = g2) + 
  wrap_elements(full = shared_legend) +
  plot_layout(ncol = 3, widths = c(1, 1, 0.3)) + 
  plot_annotation(title = glue('{state}')) + 
  theme(plot.title = element_text(
    family = "Arial",      # Font family
    face = "bold",         # Font style ("plain", "italic", "bold", "bold.italic")
    size = 12,             # Font size in points
    hjust = 0.01,           # Centering: 0 = Left, 0.5 = Centered, 1 = Right
    color = "grey20",      # Font color
    margin = margin(b = 10)# Spacing/margin underneath the title (bottom margin)
  ))

# 3. Arrange using grid.arrange, which natively respects ggtern grobs
pdf(glue('figs/{state}.pdf'), width = 8, height = 4)
print(p_final)
dev.off()


