data_env <- read_rds('output/glioma_condition_data_env.rds')
list2env(as.list(data_env), envir = globalenv())
res.ARGO <- read_rds('output/WT_vs_mutant_condition_results.rds')
funGP_fits <- res.ARGO$resObjects

funGP_res <- funGP_fits %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation)) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 

#============================================================
#============================================================
v_thresh <- 1
h_thresh <- -log(0.05)
p <- funGP_res %>% mutate(status = case_when(
  neglogp > h_thresh & deviation > v_thresh ~ 'Upreg',
  neglogp > h_thresh & deviation < -v_thresh ~ 'Downreg',
  TRUE ~ 'NS'
)) %>% arrange(status == 'Not Significant') %>%
  ggplot(aes(x = deviation, y = neglogp)) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status, alpha = status), size = 2.5, shape = 20) +
  scale_color_manual( name = "",
                      values = c(
                        "Upreg" = "#e41a1c",   # Professional Red
                        "Downreg" = "#377eb8",  # Professional Blue
                        "NA" = "gray60"     # Muted Gray
                      )
  ) + scale_alpha_manual( name = NULL, 
                          guide = "none",
                          values = c(
                            "Upreg" = .5,   # Professional Red
                            "Downreg" = .5,  # Professional Blue
                            "NS" = .1     # Muted Gray
                          )) + 
  geom_text_repel(
    data = . %>% arrange(desc(stats), desc(deviation)), #%>%
      # {bind_rows(slice_head(., n=5), slice_tail(., n=))},
    aes(label = gene),
    size = 8 / .pt,
    family = 'Arial',
    color = "black",
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "gray", 
    max.overlaps = 20
  ) +
  # 2. Horizontal significance threshold line: -log10(0.05) ~ 1.301
  geom_hline(yintercept = h_thresh, color = "red", linetype = "dashed", linewidth = .75) +
  # 3. Vertical effect size threshold lines
  geom_vline(xintercept = c(-v_thresh, v_thresh), color = "red", linetype = "dashed", linewidth = .75) +
  # 4. Axis Labels
  labs(
    x = "max deviation (LFC)",
    y = expression(-log10~" analytic pval")
  ) + scale_y_log10() + 
  # 5. Clean, professional theme minimalization
  theme_bw(base_family = 'Arial', 
           base_size = 12) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size=8),
  ) + theme(
    legend.position = "right", 
    # legend.title = element_text(size = 12, face = "bold"), # Title font size
    legend.text = element_text(size = 8),
    legend.box = "horizontal",
    legend.margin = margin(l = -5, r = 0, b = 0, t = 0), # Pulls legend closer to plot
    # legend.text = element_text(size = 10),              # Reduces font size if needed
    legend.key.size = unit(3, "mm")                 # Shrinks the color boxes
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1))
  ) + scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  # scale_y_continuous(
  #   expand = expansion(mult = c(0.025, 0.15))) +
  coord_cartesian(clip = "off")

print(p)
ggsave(plot = p, filename = glue('figs/{KERNEL}_{CONTRAST}_abundance_volcano.pdf'), 
       width = 4, height = 3, units = 'in')

#============================================================
#============================================================

# visualize abundances in 3D
all(metadata$grid_key == rownames(counts))
bin_meta <- cbind(metadata, counts)

#============================================================
#============================================================
library(scatterplot3d)

col = ifelse(bin_meta$condition == "IDHmutant", "#1f77b4", "#d62728")
col_transparent <- adjustcolor(col, alpha.f = 0.5)
psz <- scales::rescale(bin_meta$cMESInf, to = c(0.5, 3))

dxf <- bin_meta %>% select(bin, CT, CTmvp, CTpan, LE) %>% distinct() %>% arrange(bin)
pan <- dxf$CTpan
le <- dxf$LE
mvp <- dxf$LE

## need to highlight ROI too!
stats <- funGP_fits$cMESInf$statistic
band <- funGP_fits$cMESInf$band
idx <- which((stats > band) | (stats < -band) )

plot(le, stats, ylim = range(c(stats, band, -band)), col = 'blue', pch=20, cex=1.5)
points(le, band, col = 'red', pch = 20)
points(le, -band, col = 'red', pch = 20)


s3d <- scatterplot3d(
  x = bin_meta$CTmvp, 
  y = bin_meta$CTpan, 
  z = bin_meta$LE,
  pch = 16,                  # Solid circle markers
  color = col_transparent,          # Color by condition
  cex.symbols = psz, # Scale marker size by abundance
  main = "cMesInf",
  xlab = "CTmvp",
  ylab = "CTpan",
  zlab = "LE",
  angle = 45,                # Adjust this angle (in degrees) to rotate the x-y plane
  grid = FALSE,
  box = TRUE
)

legend("topright",
       legend = c("Mutant", "WT"),
       col = c("#1f77b4", "#d62728"),
       pch = 16,
       pt.cex = 1.25,         # Fixed size for the legend points
       inset = 0.05,
       bty = "n")          # No box around the legend

