mean_abx <- rowsum(prop, bin_meta$condition_key) %>% 
  as.data.frame() %>% rownames_to_column('condition_key') 

mean_abx <- pivot_longer(mean_abx, cols = -condition_key, 
                         names_to = 'cell_type', 
                         values_to = 'proportion') %>% 
  left_join(bin_meta %>% select(condition, coord_x, condition_key),
            by = 'condition_key', multiple = 'first')

## only plot the olig cell types
mean_abx <- filter(mean_abx, cell_type %in% c("Olig-1", "Olig-2", "Olig-3"))

pal43 <- colorRampPalette(RColorBrewer::brewer.pal(6, "Set2"))(3)
p <- ggplot(mean_abx, aes(x = coord_x, 
                          y = proportion,  # Or n_cells
                          color = cell_type, 
                          group = cell_type)) +
  geom_line(linewidth = 1, alpha = 0.75, linetype=2) +
  facet_wrap(~ interaction(condition), 
             scales = "free_y", ncol=1) +
  labs(
    # title = "Mean Cell Type Proportions vs. Isodepth Bins by Age",
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
        axis.text = element_text(size = 6, face = "plain"),
        
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
  

ggsave(plot = p, filename = glue('figs/Olig_fractions.pdf'), 
       width = 4, height =4, units = 'in')
