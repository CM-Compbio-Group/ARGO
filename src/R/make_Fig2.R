funGP_fits <- res.ARGO$resObjects

res <- filter(funGP_res, abs(net_poe) > .25 & abs(net_poe) < .6) %>% 
  filter(fdr < .05) %>% arrange(desc(neglogp)) %>% 
  slice_head(., n=30)
  # {bind_rows(slice_head(., n = 20), slice_tail(., n = 10))}

stat_lfc <- compute_stat_continuous(funGP_fits, N)
col_fun = colorRamp2(
  breaks = c(-1, 0, 1),                        # Min, Center, Max limits
  colors = c("#377eb8", "lightgrey", "#e41a1c")    # Blue to Red gradient
)

pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_partial_heatmap.pdf'), width = 4, height =4)
Heatmap(stat_lfc[res$gene, ], cluster_columns = F, name = 'LFC',
        show_row_names = T,
        show_column_names = F, 
        col = col_fun,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8), 
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 8, fontface = "bold"), # Title size 14 pt
          labels_gp = gpar(fontsize = 8)                    # Labels size 10 pt
        ))
dev.off()



#=====================================================
#=====================================================
mean_abx <- rowsum(counts, metadata$condition_key) %>% 
  as.data.frame() %>% rownames_to_column('condition_key') 

mean_abx <- pivot_longer(mean_abx, cols = -condition_key, 
                         names_to = 'cell_type', 
                         values_to = 'countsortion') %>% 
  left_join(metadata %>% select(condition, coord_x, condition_key),
            by = 'condition_key', multiple = 'first')

## only plot the olig cell types
mean_abx <- filter(mean_abx, cell_type %in% c("Olig-1", "Olig-2", "Olig-3"))

pal43 <- colorRampPalette(RColorBrewer::brewer.pal(6, "Set2"))(3)
pal43 <- c('gray85', 'slategray', 'gray25')
pal43 <- c('lightpink', 'palevioletred', 'darkred')

p <- ggplot(mean_abx, aes(x = coord_x, 
                          y = countsortion,  # Or n_cells
                          color = cell_type, 
                          group = cell_type)) +
  geom_line(linewidth = 1, alpha = 0.75, linetype=2) +
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
       width = 4, height =4, units = 'in')
