library(ggtern)
library(dplyr)
library(tidyr)

# 2. Helper function to generate gated Differential Ternary Plot for any program
plot_differential_niche <- function(program_name, roi_bins, delta_thresh = 0.25, active_thresh = 0.25) {
  
  # Build paired table per bin: WT vs Mutant
  df_diff <- bmf %>% 
    group_by(condition_key) %>% 
    mutate(val = mean(.data[[program_name]], na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(bin, condition, simplex_CTpan, simplex_CTmvp, simplex_LE, val) %>%
    pivot_wider(names_from = condition, values_from = val) %>%
    mutate(
      diff_WT_vs_Mut = GBM - IDHmutant,
      max_activity   = pmax(GBM, IDHmutant),
      inROI = ifelse(
        (bin %in% roi_bins) & (abs(diff_WT_vs_Mut) > delta_thresh) & (max_activity > active_thresh),
        "T", "F"
      )
    ) %>%
    arrange(inROI == "T")
  
  ggtern(data = df_diff, aes(x = simplex_CTpan, y = simplex_CTmvp, z = simplex_LE)) +
    geom_point(
      aes(fill = diff_WT_vs_Mut, size = inROI, alpha = inROI, stroke = inROI),
      shape = 21, color = "black"
    ) +
    # Hide major/minor grid lines and stray L/R/T indicator arrows
    theme_nogrid() +
    theme_hidearrows() +
    scale_fill_gradient2(
      low = "#2166ac",      # Blue = IDHmutant
      mid = "#f7f7f7",      # Grey-white = No difference
      high = "#b2182b",     # Crimson = IDHwt
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = bquote(atop(Delta ~ .(program_name), "(WT - Mut)")),
      guide = guide_colorbar(
        title.position = "top",
        barheight = unit(0.8, "in"),
        barwidth = unit(0.12, "in"),
        frame.colour = "black"
      )
    ) +
    scale_size_manual(values = c("F" = 1.2, "T" = 2), name = "Active ROI") +
    scale_alpha_manual(values = c("F" = 0.45, "T" = 0.95), guide = "none") +
    scale_discrete_manual(aesthetics = "stroke", values = c("F" = 0.3, "T" = 0.8), guide = "none") +
    labs(
      title = program_name,
      x = "CTpan", y = "CTmvp", z = "LE"
    ) +
    theme(
      # Universal text family and size
      text = element_text(family = "Arial", size = 8),
      
      # Ternary axes & tick marks
      tern.axis.line = element_line(color = "black", linewidth = 0.5),
      tern.axis.text = element_text(family = "Arial", size = 8),
      tern.axis.title = element_text(family = "Arial", size = 8),
      tern.panel.background = element_blank(),
      
      # Plot title
      plot.title = element_text(family = "Arial", size = 12, face = "bold", 
                                hjust = 0, margin = margin(b = 2)),
      
      # Legend titles and labels
      legend.title = element_text(family = "Arial", size = 8),
      legend.text = element_text(family = "Arial", size = 8),
      
      # Tight spacing
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.spacing = unit(3, "pt")
    )
}
#============================================================
#============================================================
## visualize posterior fitted counts in significant bins
bmf <- left_join(full_grid, bin_meta %>% distinct(sample_id, condition), by = 'sample_id') %>%
  mutate(condition_key =  paste(condition, bin, sep = "___")) %>%
  mutate(condition_key = factor(condition_key, 
                                levels = gtools::mixedsort(unique(condition_key)))
         )

data_env$features_to_fit <- c('Vasc', 'cMESHyp', 'cUndiff', 'cProlifMetab')

data_env$features_to_fit <- c('cMESVasc', 'cOPC2', 'Synaptic')

data_env$kernel_params <- res.ARGO$params
data_env$GET_POSTERIOR <- TRUE
y_fit <- fit_ARGO(data_env)
y_fit <- do.call('cbind', y_fit)
colnames(y_fit) <- data_env$features_to_fit
bmf <- bind_cols(bmf, y_fit)


roi_index <- lapply(funGP_fits, function(kxx) which(abs(kxx$statistic) > kxx$band))
sapply(roi_index, length)

roi_index$Vasc

# brain <- read_rds('~/Desktop/brain-seurat-objects/visium.rds')
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

bmf <- left_join(bmf, simplex_df, by = "bin")
for(celltype in data_env$features_to_fit){
  message(celltype)
  p <- plot_differential_niche(program_name = celltype, roi_bins = roi_index[[celltype]])
  ggsave(plot = p, filename = glue('figs/{celltype}_ternary_delta.pdf'), width = 3.5, height=3.5)
}


#==========================================================
#==========================================================

bin_meta %>% group_by(condition_key) %>% 
  mutate(val = mean(.data[[program_name]]), na.rm=T) %>%
  ungroup() %>%
  distinct(bin, condition, val) %>% filter(bin %in% roi_bins) %>% 
  pivot_wider(names_from = condition, values_from = val)

bin_meta %>% group_by(condition_key) %>% 
  mutate(val = mean(.data[[program_name]]), na.rm=T) %>% pull(val) %>% plot()


# ggsave(plot=p, filename = 'figs/delta_abundance.pdf', width = 8, height = 8)
