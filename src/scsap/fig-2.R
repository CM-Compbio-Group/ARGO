library(Seurat)
dxlist <- read_rds('~/Desktop/slide_list.rds')
dxlist

## visualize the correspondence b/w isodepth and tissue domains in these slides

title_text <- c('4wk (juvenile)', '24wk (adult)', '90wk (old)')
my_plots <- Map(function(dx, title_text) {
  FeaturePlot(dx, features = 'isodepth',
              # group.by = 'tissue',
          reduction = 'X_spatial_coords',
          alpha = .75,
          # pt.size = .75
          ) +
    scale_color_gradientn(
      colors = c("lightgray",  "red2", "darkred"),
      limits = c(0,1),
    ) +
  # DimPlot(dx, group.by = 'tissue', reduction = 'X_spatial_coords', alpha = .75)+
    ggtitle(title_text) + # <-- Your custom title goes here dynamically
    xlab('') + ylab('') + labs(color = 'Isodepth') + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line = element_blank())
}, dxlist[c(1, 3, 5)], title_text)



p <- wrap_plots(my_plots, ncol = 3) + 
  plot_layout(guides = "collect") 


ggsave(plot = p, filename = 'figs/domain_isodepth.pdf', width = 12, height=4)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

## show mapping between tissue domains and isodepth across all slides
# distribution of isodepth values across tissue domains by slide
obs <- read_rds('data/temp/ad_obs.rds')
obs <- obs %>% mutate(tissue = fct_relevel(tissue, 'striatum', 'brain ventricle'))

p <- ggplot(obs, aes(x = factor(tissue), y = binned_isodepth, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + # Custom color palette
  theme_classic() +
  labs(x = "Tissue domain", y = "Isodepth") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=12),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size=14)
    ) + labs(fill = 'Tissue Annotation')

ggsave(plot=p, filename = 'figs/isodepth_vs_tissue_boxplots.pdf', width = 8, height = 6)  


#----------------------------------------------------------------------
#----------------------------------------------------------------------

