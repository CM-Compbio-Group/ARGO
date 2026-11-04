source('src/R/summarize_argo_excursions.R')

N_BINS <- 65
data_env <- read_rds(glue('output/glioma_kmeans_bysample_N_{N_BINS}_data_env.rds'))
# list2env(as.list(data_env), envir = globalenv())
res.ARGO <- read_rds(glue('output/glioma_kmeans_bysample_N_{N_BINS}_results.rds'))

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
#============================================================
h_thresh <- log2(1.5); v_thresh <- -log10(0.05)
p <- funGP_res %>% 
  mutate(status = case_when(
    neglogp > v_thresh & directed_deviation > h_thresh ~ 'Upreg',
    neglogp > v_thresh & directed_deviation < -h_thresh ~ 'Downreg',
    TRUE ~ 'NS'
  )) %>% 
  arrange(status == 'NS') %>%
  ggplot(aes(x = directed_deviation, y = neglogp)) +
  geom_point(aes(color = status, alpha = status, size = polarity), shape = 20) +
  scale_color_manual(
    name = bquote(Delta * " abundance"),
    values = c(
      "Upreg" = "#e41a1c",
      "Downreg" = "#377eb8",
      "NS" = "gray70"
    )
  ) + 
  scale_alpha_manual(
    name = "", 
    guide = "none",
    values = c("Upreg" = 0.5, "Downreg" = 0.5, "NS" = 0.1)
  ) + 
  scale_size_continuous(
    name = 'Polarity',
    range = c(.25, 4),
    limits = c(0.5, 1),
    trans = "identity"
  ) + 
  geom_text_repel(
    data = . %>% arrange(desc(stats), desc(directed_deviation)),
    aes(label = feature),
    size = 6 / .pt,
    family = 'Arial',
    color = "black",
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "gray", 
    max.overlaps = 20
  ) +
  geom_hline(yintercept = v_thresh, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = c(-h_thresh, h_thresh), color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(
    x = "Max directed deviation",
    y = expression(-log[10]~" analytic p-value")
  ) + 
  theme_bw(base_family = 'Arial', base_size = 8) +
  theme(
    text = element_text(family = 'Arial', size = 8),
    axis.title = element_text(family = 'Arial', size = 8),
    axis.text = element_text(family = 'Arial', size = 8),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # Legend layout & compact vertical spacing
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(family = 'Arial', size = 8),
    legend.text = element_text(family = 'Arial', size = 8),
    
    # Compress key height and vertical spacing between rows
    legend.key.width   = unit(3.0, "mm"),
    legend.key.height  = unit(2.5, "mm"),
    legend.key.spacing.y = unit(0, "pt"),     # ggplot2 >= 3.5.0
    legend.spacing.y   = unit(2, "pt"),       # Gap between color & size boxes
    
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.spacing = unit(3, "pt"),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  ) + 
  guides(
    color = guide_legend(
      order = 1,
      keyheight = unit(2, "mm"),
      override.aes = list(size = 3.5, alpha = 1)
    ), 
    size = guide_legend(
      order = 2,
      keyheight = unit(2, "mm")             # Locks row height for large points
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  coord_cartesian(clip = "off")

print(p)

ggsave(plot = p, filename = glue('figs/{data_env$KERNEL}_{data_env$CONTRAST}_N_{N_BINS}_abundance_volcano.pdf'),
       width = 4, height = 3, units = 'in')

#============================================================
#============================================================
# # visualize abundances in 3D
# all(bin_meta$grid_key == rownames(data_env$counts))
# bin_meta <- cbind(bin_meta, counts)
# #============================================================
# dxf <- bin_meta %>% select(bin, CT, CTmvp, CTpan, LE) %>% distinct() %>% arrange(bin)
# pan <- dxf$CTpan
# le <- dxf$LE
# mvp <- dxf$LE
# 
# ## need to highlight ROI too!
# stats <- funGP_fits$cMESInf$statistic
# band <- funGP_fits$cMESInf$band
# idx <- which((stats > band) | (stats < -band) )
# plot(le, stats, ylim = range(c(stats, band, -band)), col = 'blue', pch=20, cex=1.5)
# points(le, band, col = 'red', pch = 20)
# points(le, -band, col = 'red', pch = 20)
# #============================================================
# library(scatterplot3d)
# 
# col = ifelse(bin_meta$condition == "IDHmutant", "#1f77b4", "#d62728")
# col_transparent <- adjustcolor(col, alpha.f = 0.5)
# psz <- scales::rescale(bin_meta$cMESInf, to = c(0.5, 3))
# 
# s3d <- scatterplot3d(
#   x = bin_meta$CTmvp, 
#   y = bin_meta$CTpan, 
#   z = bin_meta$LE,
#   pch = 16,                  # Solid circle markers
#   color = col_transparent,          # Color by condition
#   cex.symbols = psz, # Scale marker size by abundance
#   main = "cMesInf",
#   xlab = "CTmvp",
#   ylab = "CTpan",
#   zlab = "LE",
#   angle = 45,                # Adjust this angle (in degrees) to rotate the x-y plane
#   grid = FALSE,
#   box = TRUE
# )
# 
# legend("topright",
#        legend = c("Mutant", "WT"),
#        col = c("#1f77b4", "#d62728"),
#        pch = 16,
#        pt.cex = 1.25,         # Fixed size for the legend points
#        inset = 0.05,
#        bty = "n")          # No box around the legend

