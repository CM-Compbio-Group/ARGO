source('src/argo/00-load-packages.R')

data_env <- read_rds('output/abundance_data_env.rds')

data_env$CONTRAST <- 'old_vs_adult'
data_env$CONTRAST <- 'adult_vs_juvenile'

list2env(as.list(data_env), envir = globalenv())

res.ARGO <- read_rds(glue('output/{data_env$CONTRAST}_diff_abundance_results.rds'))
funGP_fits <- res.ARGO$resObjects

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
## convert to mean scaled counts in each conditionxbin 
n_counts <- metadata %>% dplyr::count(condition_key)
mean_abx <- rowsum(counts, metadata$condition_key) 
## divide by frequency of condition key in metadata
all(rownames(mean_abx) == n_counts$condition_key)
mean_abx <- mean_abx / n_counts$n 

mean_abx <- mean_abx %>% 
  as.data.frame() %>% 
  rownames_to_column('condition_key') 

mean_abx <- pivot_longer(mean_abx, cols = -condition_key, 
                         names_to = 'cell_type', 
                         values_to = 'counts') %>% 
  left_join(metadata %>% select(condition, coord_x, condition_key),
            by = 'condition_key', multiple = 'first')

## only plot the olig cell types
mean_abx <- filter(mean_abx, cell_type %in% c("Olig-1", "Olig-2", "Olig-3"))

pal43 <- colorRampPalette(RColorBrewer::brewer.pal(6, "Set2"))(3)
pal43 <- c('gray85', 'slategray', 'gray25')
pal43 <- c('lightpink', 'palevioletred', 'darkred')

p <- ggplot(mean_abx, aes(x = coord_x, 
                          y = counts,  # Or n_cells
                          color = cell_type, 
                          group = cell_type)) +
  geom_line(linewidth = .75, alpha = 0.75, linetype=2) +
  facet_wrap(~ interaction(condition), 
             scales = "free_y", ncol=1) +
  labs(
    # title = "Mean Cell Type countsortions vs. Isodepth Bins by Age",
    x = "Isodepth", y = "Scaled Abundance",
    color = "Cell Type") +
  scale_color_manual(values = pal43) +
  theme_classic(base_family = 'Arial', 
                base_size = 8) +
  # guides(color = guide_legend(ncol = 3)) +  # Set number of columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom", 
        plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        
        strip.text = element_text(size = 8, face = "plain"), 
        # strip.background = element_blank(),           # Optional: Cleans up facet header boxes
        
        # Legend text and title formatting
        legend.title = element_text(size = 8, face = "plain"),
        legend.text = element_text(size = 8, face = "plain"), 
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0), # Pulls legend closer to plot
        legend.key.size = unit(3, "mm")                 # Shrinks the color boxes
  ) + scale_y_continuous(
    expand = expansion(mult = c(0.025, 0.15))) +
  coord_cartesian(clip = "off") 

print(p)
ggsave(plot = p, filename = glue('figs/Olig_fractions.pdf'), 
       width = 2.5, height =3, units = 'in')
