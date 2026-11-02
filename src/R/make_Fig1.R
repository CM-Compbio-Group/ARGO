source('src/argo/00-load-packages.R')
source('src/R/summarize_argo_excursions.R')

# data_env <- read_rds('output/expression_data_env.rds')
# list2env(as.list(data_env), envir = globalenv())
res.ARGO <- read_rds(glue('output/old_vs_adult_N30_wo_donor_diff_testing.rds'))

funGP_fits <- res.ARGO$resObjects
funGP_summaries <- lapply(funGP_fits, summarize_argo_excursions)
funGP_res <- do.call('rbind', funGP_summaries) %>% 
  mutate(feature = names(funGP_fits), 
         fdr = p.adjust(exp(-neglogp), method = 'BH'),
         stats = neglogp * net_sign) %>%
  relocate(feature, .before = everything()) %>% 
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())
print(funGP_res)

#============================================================
## how many significant genes
filter(funGP_res, fdr < .05 & abs(directed_deviation) > 1) %>% 
  dplyr::count(net_sign)
#============================================================
#============================================================
h_thresh <- log2(2); v_thresh <- -log10(0.05)
n_label <- 10
p <- funGP_res %>%
  mutate(status = case_when(
    -log10(fdr) > v_thresh & directed_deviation > h_thresh ~ 'Upreg',
    -log10(fdr) > v_thresh & directed_deviation < -h_thresh ~ 'Downreg',
    TRUE ~ 'NS'
  )) %>%
  arrange(status == 'Not Significant') %>%
  ggplot(aes(x = directed_deviation, y = -log10(fdr))) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status, alpha = status, size = polarity), shape = 20) +
  scale_color_manual(
    name = bquote(Delta * " expression"),
    values = c(
      "Upreg" = "#e41a1c",   # Professional Red
      "Downreg" = "#377eb8", # Professional Blue
      "NA" = "gray70"        # Muted Gray (Fixed typo from "NA" to "NS" based on your mutate)
    )
  ) +
  scale_alpha_manual(name = "", guide = "none", values = c(
    "Upreg" = .35,
    "Downreg" = .35,
    "NS" = .1
  )) +
  scale_size_continuous(
    name = 'Polarity',
    range = c(.25, 2),
    limits = c(0.5, 1),
    trans = "identity"
  ) +
  geom_text_repel(
    data = . %>% 
      arrange(desc(stats), desc(directed_deviation)) %>% 
      {bind_rows(slice_head(., n=n_label), slice_tail(., n=n_label))},
    aes(label = feature),
    size = 8 / .pt,
    family = 'Arial',
    color = "black",
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "gray",
    max.overlaps = 20
  ) +
  # 2. Horizontal significance threshold line: -log10(0.05) ~ 1.301
  geom_hline(yintercept = v_thresh, color = "red", linetype = "dashed", linewidth = .75) +
  # 3. Vertical effect size threshold lines
  geom_vline(xintercept = c(-h_thresh, h_thresh), color = "red", linetype = "dashed", linewidth = .75) +
  # 4. Axis Labels
  labs(
    x = "Max directed deviation",
    y = expression(-log10~" FDR")
  ) +
  # 5. Clean, professional theme minimalization
  theme_bw(base_family = 'Arial', base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    
    # --- FIXES FOR MARGINS AND SPACING ---
    
    # 1. Reduces the outer margin of the entire plot (Top, Right, Bottom, Left)
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"), 
    
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.box = "vertical",
    
    # 2. Reduces padding inside the legend box itself
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
    
    # 3. Pulls the entire legend box closer to the plot panel (negative left margin)
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -5, unit = "pt"), 
    
    legend.key.size = unit(3, "mm")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 5, alpha = 1)),
    size = guide_legend(order = 2)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  coord_cartesian(clip = "off")

print(p)

ggsave(plot = p, filename = glue('figs/{data_env$KERNEL}_{data_env$CONTRAST}_volcano.pdf'),
       width = 3.5, height = 2.5, units = 'in')


#============================================================
#============================================================

features_to_plot <- funGP_res %>% arrange(desc(stats)) %>% 
  {bind_rows(slice_head(., n = 2), 
             slice_tail(., n = 2))} %>% pull(feature)


data_env$kernel_params <- res.ARGO$params
data_env$features_to_fit <- features_to_plot
data_env$PLOT_POSTERIOR <- T
data_env$PLOT_CONF_BANDS <- T
data_env$RETURN_FITS <- F

fit_ARGO(data_env)

#============================================================
#============================================================
require(readxl)
require(UpSetR)
library(PRROC)

res.list <- read_rds(glue('output/{data_env$CONTRAST}_DGE_results_list.rds'))
res.list$ARGO <- funGP_res %>% dplyr::rename(gene = feature)

res.sig <- lapply(res.list, function(res) filter(res, fdr < 0.05) %>% pull(gene))
names(res.sig) <- sprintf("  %-10s", names(res.sig))
sapply(res.sig, length)

pdf('figs/upsetPlot_FDR_05.pdf', width = 3, height = 2.5)
theme_set(theme_gray(base_family = "Arial", 
                     base_size = 8))
# 2. Render plot with tuned text scales and matrix ratio
print(
  upset(
    fromList(res.sig), 
    order.by = "freq", 
    nsets = length(res.sig),
    nintersects = 8, 
    shade.alpha = 0.5, 
    point.size = 2, 
    line.size = 0.75,
    
    # Allocates clean vertical room between top bars and bottom matrix
    mb.ratio = c(0.65, 0.35),
    
    # text.scale order:
    # 1: Intersection size bar numbers (top)
    # 2: Intersection size axis title (Y-axis)
    # 3: Set size bar labels (numbers on bars)
    # 4: Set size axis title (X-axis)
    # 5: Set names (rows on grid) -> slightly smaller helps prevent clipping
    # 6: Numbers above intersection bars
    text.scale = c(0.75, 0.7, 0.45, 0.6, 0.65, 0.75)
  ), 
  newpage = FALSE
)

dev.off()



funGP_res <- funGP_res %>% mutate(net_frac_bin = (n_pos_bins + n_neg_bins) / 30)
res <- filter(funGP_res, net_frac_bin > .2 & net_frac_bin < .65) %>% 
  filter(fdr < .05) %>% arrange(desc(neglogp)) %>% 
  slice_head(., n=40)
# {bind_rows(slice_head(., n = 20), slice_tail(., n = 10))}

stat_lfc <- compute_stat_continuous(funGP_fits, N)
col_fun = colorRamp2(
  breaks = c(-1, 0, 1),                        # Min, Center, Max limits
  colors = c("#377eb8", "lightgrey", "#e41a1c")    # Blue to Red gradient
)

pdf(glue('figs/{data_env$KERNEL}_{data_env$CONTRAST}_FC_partial_heatmap.pdf'), width = 4.5, height =4)
Heatmap(stat_lfc[res$feature, ], cluster_columns = F, 
        name = "delta",
        show_row_names = T,
        show_column_names = F, 
        col = col_fun,
        row_names_gp = gpar(fontsize=6),
        column_names_gp = gpar(fontsize=8), 
        heatmap_legend_param = list(
          title = expression(paste(Delta, " expression")), # Formatted label placed here
          title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
        ))
dev.off()
