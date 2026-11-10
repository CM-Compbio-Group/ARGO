source('src/argo/00-load-packages.R')
source('src/R/summarize_argo_excursions.R')

data_env <- read_rds('output/abundance_data_env.rds')
# list2env(as.list(data_env), envir = globalenv())

data_env$CONTRAST <- 'old_vs_adult'
data_env$CONTRAST <- 'adult_vs_juvenile'

res.ARGO <- read_rds(glue('output/{data_env$CONTRAST}_diff_abundance_results.rds'))

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
## i want confidence bands
data_env$CONTRAST <-  'adult_vs_juvenile'
data_env$CONTRAST_VECTOR <- c(-1, 1,0)
data_env$features_to_fit <- c('Olig-1','Olig-2','Olig-3')
data_env$kernel_params <- res.ARGO$params
data_env$GET_POSTERIOR <- FALSE
data_env$PLOT_POSTERIOR <- FALSE
data_env$PLOT_CONF_BANDS <- TRUE
data_env$RETURN_FITS <- FALSE

ci_plots <- fit_ARGO(data_env)
for(fx in data_env$features_to_fit){
  message(fx)
  ggsave(plot = ci_plots[[fx]], filename = glue('figs/{data_env$CONTRAST}_{fx}_ci.pdf'), 
         width=2.5, height=1.5, units='in')
}

#============================================================
#============================================================
h_thresh <- log2(1.5); v_thresh <- -log10(0.05)
p <- funGP_res %>% mutate(status = case_when(
  -log10(fdr) > v_thresh & directed_deviation > h_thresh ~ 'Upreg',
  -log10(fdr) > v_thresh & directed_deviation < -h_thresh ~ 'Downreg',
  TRUE ~ 'NS'
)) %>% arrange(status == 'Not Significant') %>%
  ggplot(aes(x = directed_deviation, y = -log10(fdr))) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status, alpha = status, size = polarity),  shape = 20) +
  scale_color_manual( name = bquote(Delta * " abundance"),
                      values = c(
                        "Upreg" = "#e41a1c",   # Professional Red
                        "Downreg" = "#377eb8",  # Professional Blue
                        "NA" = "gray70"     # Muted Gray
                      )
  ) + scale_alpha_manual(name = "", 
                         guide = "none",
                         values = c(
                           "Upreg" = .5,   # Professional Red
                           "Downreg" = .5,  # Professional Blue
                           "NS" = .1     # Muted Gray
                         )) + 
  scale_size_continuous(
    name = 'Polarity',
    range = c(1, 4),
    # breaks = waiver(),
    # labels = waiver(),
    limits = c(0.5, 1),
    trans = "identity"
  ) + 
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) + 
  geom_text_repel(
    data = . %>% arrange(desc(stats), desc(directed_deviation)) %>%
    {bind_rows(slice_head(., n=5), slice_tail(., n=5))},
    aes(label = feature),
    size = 8 / .pt,
    family = 'Arial',
    color = "black",
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "gray", 
    max.overlaps = 30
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
  # scale_y_log10() + 
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
    legend.box = "vertical",
    legend.margin = margin(t = 7.5, r = 2.5, b = 5, l = 2.5, unit = "mm"), # Pulls legend closer to plot
    # legend.text = element_text(size = 10),              # Reduces font size if needed
    legend.key.size = unit(3, "mm")                 # Shrinks the color boxes
  ) + 
  guides(
    color = guide_legend(order = 1,override.aes = list(size = 5, alpha = 1)), 
    size  = guide_legend(order = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  # scale_y_continuous(
  #   expand = expansion(mult = c(0.025, 0.15))) +
  coord_cartesian(clip = "off")

print(p)

ggsave(plot = p, filename = glue('figs/{data_env$KERNEL}_{data_env$CONTRAST}_abundance_volcano.pdf'), 
       width = 4, height = 2.5, units = 'in')


exc.metrics <- do.call('rbind', lapply(funGP_fits, compute_exceedence_metrics, isodepth_))
## aggregate results from functional testing
funGP_res <- funGP_fits %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  # bind_cols(exc.metrics) %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         # net_auec = auec_up - auec_down,
         # net_poe = poe_up - poe_down
  ) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 

funGP_res

#=====================================================
#=====================================================

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
    legend.position = "none", 
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

print(p)
ggsave(plot = p, filename = glue('figs/{KERNEL}_{CONTRAST}_abundance_volcano.pdf'), 
       width = 2.5, height = 2.5, units = 'in')


#=====================================================
#=====================================================

