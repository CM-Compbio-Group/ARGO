source('src/argo/00-load-packages.R')

KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'

data_env <- read_rds('output/expression_data_env.rds')

data_env$cov_donor <- data_env$donor_identity
data_env$PLOT_POSTERIOR <- TRUE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$RETURN_FITS <- FALSE
data_env$features_to_fit <- features_to_plot
data_env$kernel_params <- res.ARGO$params # if my parameters are already fit, I don't want to fit likelihood again

fit_ARGO(data_env)


features_to_plot <- funGP.df %>% arrange(desc(stats)) %>% 
  {bind_rows(slice_head(., n = 2), 
             slice_tail(., n = 2))} %>% pull(gene)

res.ARGO <- read_rds(glue('output/{CONTRAST}_N30_diff_testing_results.rds'))
funGP.df <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())


v_thresh <- 1
h_thresh <- -log(0.05)

p <- funGP_res %>% mutate(status = case_when(
  neglogp > h_thresh & deviation > v_thresh ~ 'Upreg',
  neglogp > h_thresh & deviation < -v_thresh ~ 'Downreg',
  TRUE ~ 'NS'
)) %>% arrange(status == 'Not Significant') %>%
  ggplot(aes(x = deviation, y = -log10(fdr))) +
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
    data = . %>% arrange(desc(stats), desc(deviation)) %>% {bind_rows(slice_head(., n=5), slice_tail(., n=5))},
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
  geom_vline(xintercept = c(-v_thresh, v_thresh), color = "red", linetype = "dashed", size = .75) +
  # 4. Axis Labels
  labs(
    x = "max deviation (LFC)",
    y = expression(-log10~" FDR")
  ) +
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
  scale_y_continuous(
    expand = expansion(mult = c(0.025, 0.15))) +
  coord_cartesian(clip = "off")

ggsave(plot = p, filename = glue('figs/{data_env$KERNEL}_{data_env$CONTRAST}_abundance_volcano.pdf'), 
      width = 4, height = 3, units = 'in')


require(readxl)
require(UpSetR)
library(PRROC)

res.list <- read_rds(glue('output/{CONTRAST}_DGE_results_list.rds'))
res.sig <- lapply(res.list, function(res) filter(res, fdr < 0.05) %>% pull(gene))
sapply(res.sig, length)

pdf('figs/upsetPlot_FDR_05.pdf', width = 4, height = 4)
theme_set(theme_gray(base_family = "Arial", 
                     base_size = 8))
print(upset(fromList(res.sig), 
            order.by = "freq", 
            nsets = length(res.sig),
            nintersects = 10, shade.alpha = .5, point.size = 2, line.size = .75,
            text.scale = c(
              .75, # 1. Intersection size bar numbers (top)
              .5, # 2. Intersection size axis title (Y-axis label)
              .25, # 3. Set size bar labels (left)
              .5, # 4. Set size axis title (X-axis label)
              .75, # 5. Set names (rows on the grid matrix)
              .75  # 6. Numbers above the dots (if turned on)
            )), 
      newpage= F)
dev.off()


