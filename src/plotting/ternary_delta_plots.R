library(ggtern)
library(dplyr)
library(tidyr)

ivy_tern <- c('CTpan', 'CTmvp', 'LE')

# 1. Evenly spread bins across the triangle face using marginal quantiles (percent_rank)
raw_bin_coords <- brain@meta.data %>%
  mutate(bin = brain$bins) %>%
  group_by(bin) %>%
  summarise(across(all_of(ivy_tern), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

tern_mat <- raw_bin_coords %>%
  select(all_of(ivy_tern)) %>%
  as.matrix()

# Quantile scaling + small offset (0.05) distributes points across the interior face
tern_quantiles <- apply(tern_mat, 2, function(x) {
  (dplyr::percent_rank(x) + 0.05)^1.2
})

simplex_coords <- sweep(tern_quantiles, 1, rowSums(tern_quantiles), '/')
colnames(simplex_coords) <- paste('simplex', colnames(simplex_coords), sep = '_')
simplex_df <- bind_cols(bin = raw_bin_coords$bin, as.data.frame(simplex_coords))

bin_meta_clean <- bin_meta %>%
  select(-any_of(c("simplex_CT", "simplex_CTmvp", "simplex_CTpan", "simplex_LE"))) %>%
  left_join(simplex_df, by = "bin")

# 2. Helper function to generate gated Differential Ternary Plot for any program
plot_differential_niche <- function(program_name, roi_bins, delta_thresh = 0.5, active_thresh = 0.25) {
  
  # Build paired table per bin: WT vs Mutant
  df_diff <- bin_meta_clean %>%
    select(bin, condition, simplex_CTpan, simplex_CTmvp, simplex_LE, val = all_of(program_name)) %>%
    pivot_wider(names_from = condition, values_from = val) %>%
    mutate(
      # Signed difference: Positive = Higher in IDHwt (GBM), Negative = Higher in IDHmutant
      diff_WT_vs_Mut = GBM - IDHmutant,
      max_activity   = pmax(GBM, IDHmutant),
      # Gated ROI: statistically in ROI + meaningful difference + active above baseline
      inROI = ifelse(
        (bin %in% roi_bins) & (abs(diff_WT_vs_Mut) > delta_thresh) & (max_activity > active_thresh),
        "T", "F"
      )
    ) %>%
    # Plot non-ROI points first so ROI points sit cleanly on top
    arrange(inROI == "T")
  
  ggtern(data = df_diff, aes(x = simplex_CTpan, y = simplex_CTmvp, z = simplex_LE)) +
    geom_point(
      aes(fill = diff_WT_vs_Mut, size = inROI, alpha = inROI, stroke = ifelse(inROI == "T", 0.6, 0.1)),
      shape = 21, color = "black"
    ) +
    theme_nogrid_major() +
    scale_fill_gradient2(
      low = "#2166ac",      # Deep Blue = IDHmutant enriched
      mid = "#f7f7f7",      # Soft Grey-White = No difference
      high = "#b2182b",     # Deep Crimson/Magenta = IDHwt enriched
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = paste0("Δ ", program_name, "\n(WT - Mut)"),
      guide = guide_colorbar(
        title.position = "top",
        barheight = unit(1.0, "in"),
        barwidth = unit(0.1, "in"),
        frame.colour = "black"
      )
    ) +
    scale_size_manual(values = c("F" = 1.6, "T" = 3.8), name = "Active ROI") +
    scale_alpha_manual(values = c("F" = 0.45, "T" = 0.95), guide = "none") +
    theme(tern.panel.background = element_blank()) +
    labs(
      title = program_name,
      x = "CTpan", y = "CTmvp", z = "LE"
    )
}

# Example usage for your 4 programs:
p1 <- plot_differential_niche("cMESHyp",  roi_index$cMESHyp)
p2 <- plot_differential_niche("cMESVasc", roi_index$cMESVasc)
p3 <- plot_differential_niche("cOXPHOS",  roi_index$cOXPHOS)
p4 <- plot_differential_niche("cMESInf",  roi_index$cMESInf)



# p <- wrap_plots(list(p1,p2,p3,p4), ncol=2)
# ggsave(plot=p, filename = 'figs/delta_abundance.pdf', width = 8, height = 8)
