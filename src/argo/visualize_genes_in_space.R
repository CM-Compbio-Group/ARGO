#===================================================================
## read list of seurat objects with representative images
#===================================================================
library(Seurat)
library(RColorBrewer)

res.ARGO <- read_rds(glue('output/old_vs_adult_N30_diff_testing_results.rds'))
funGP_fits <- res.ARGO$resObjects

adxlist <- read_rds('~/Desktop/adxlist.rds')
dfm <- read_rds('data/temp/dfm.rds')
adxlist

## pick one slide each
adxlist <- dxlist[c(3,5)]

write_rds(adxlist, '~/Desktop/adxlist.rds')
tc <- adxlist[[1]]@reductions$X_spatial_coords@cell.embeddings 
adxlist[[1]]@reductions$X_spatial_coords@cell.embeddings <- -tc

img_title <- c('24wk', '90wk')
cols <- brewer.pal(n = length(unique(dfm$tissue)), name = "Set2")
make_spatial_plot <- function(dx, plot_title){
  DimPlot(dx, 
          group.by = 'tissue',
          # features = 'isodepth',
          reduction = 'X_spatial_coords', 
          alpha = .5, cols=cols,
          pt.size = .65, stroke.size = 0) +
    theme_void(base_family = 'Arial') + 
    # Locks the physical aspect ratio so the brain doesn't distort, and removes 5% internal padding
    coord_fixed(expand = FALSE) + 
    # scale_color_viridis_c(
    #   option = "plasma",
    #   name = "isodepth",
    #   limits = c(0, 1),
    #   guide = guide_colorbar(
    #     frame.colour = "black",
    #     ticks.colour = "black",
    #     barheight = unit(.75, "in"),
    #     barwidth = unit(0.1, "in")
    #   )
    # ) +
    ggtitle(plot_title) +
    theme(
      plot.margin = margin(0, 0, 0, 0), # Removes the individual plot margins
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 8), 
      plot.title = element_text(
        family = "Arial",
        face = "bold",
        size = 12,
        hjust = 0.5,
        color = "grey20",
        margin = margin(b = 2)
      )
    ) +
    xlab("") + ylab("") 
}

plots <- lapply(seq(2), function(i) make_spatial_plot(adxlist[[i]], img_title[i]))
p <- wrap_plots(plots, ncol = 2) + 
  plot_layout(guides = "collect") &
  theme(
    # plot.margin = margin(0.25, 0.25, 0, .25), # Removes the overarching patchwork container margins
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.key.spacing.y = unit(0.05, "cm")
  )

print(p)

ggsave(plot = p, filename = 'figs/brain_domains.pdf', width = 4, height = 2, units = 'in')

## prepare common metadata file
dfm <- do.call(rbind, lapply(adxlist, function(dx) dx@meta.data))
dfm$donor_id <- droplevels(dfm$donor_id)
dfm$age <- droplevels(dfm$age)

## need to map the spatial coordinates for these slides
obs <- read_rds('data/temp/ad_obs.rds')
obs <- obs %>% 
  mutate(bin = factor(ntile(binned_isodepth, n = N))) %>% 
  group_by(bin) %>% mutate(coord_x = mean(binned_isodepth)) %>% 
  ungroup() 
  
dfm$coord_x <- obs %>% filter(spot_id %in% rownames(dfm)) %>%
  arrange(match(spot_id, rownames(dfm))) %>% pull(coord_x)
write_rds(dfm, 'data/temp/dfm.rds')
#===================================================================
## visualize some representative genes in space
#===================================================================
# genes_to_plot <- c('Bcl11b', 'Npr3', 'Irf7', 'Ptprc', 'Foxp2')
genes_to_plot <- 'Serpina3n'
## get roi for each gene
roi <- lapply(funGP_fits[genes_to_plot], function(kxx) isodepth_[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)

filter(dfm, coord_x %in% roi$Serpina3n) %>% 
  dplyr::count(tissue) %>% 
  mutate(
    proportion = n / sum(n),         # Fraction of the total (0 to 1)
    percentage = (n / sum(n)) * 100  # Percentage of the total (0% to 100%)
  )
  

dxlist <- lapply(adxlist, function(dx) subset(dx, subset = tissue == 'olfactory region'))
lapply(dxlist, FeaturePlot, 'Serpina3n', reduction = 'X_spatial_coords')
lapply(dxlist, DimPlot, group.by = 'tissue', reduction = 'X_spatial_coords')


source('src/R/spatial_ROI_plot.R')

figs <- lapply(genes_to_plot, function(gene) {
  message(gene)
  rx <- roi[[gene]];
  lapply(dxlist, plot_roi_split, gene = gene, ROI = rx, metadata = dfm )
})

names(figs) <- genes_to_plot
plots <- lapply(genes_to_plot, function(gene) {
  message(gene)
  row_plot <- wrap_plots(figs[[gene]], ncol = 2) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{gene}')) &
    theme(
      legend.box = "horizontal",
      legend.box.just = "center",          # Aligns the centers of both bars
      legend.key.height = unit(4.5, "cm"), # Matches bar height to the rotated text height
      legend.title.position = "right",     # Ensures both titles stay on the right
      legend.title = element_text(family = "Arial", size = 8, 
                                  angle = 90, hjust = 0.5, margin = margin(l = 5)),
      legend.text = element_text(family = "Arial", size = 8),
      plot.title = element_text(
        family = "Arial",      # Font family
        face = "bold",         # Font style ("plain", "italic", "bold", "bold.italic")
        size = 12,             # Font size in points
        hjust = 0.01,           # Centering: 0 = Left, 0.5 = Centered, 1 = Right
        color = "grey20",      # Font color
        margin = margin(b = 10)# Spacing/margin underneath the title (bottom margin)
      )
    )
  wrap_elements(row_plot)
})

names(plots) <- genes_to_plot
for(gene in genes_to_plot){
  ggsave(plots[[gene]], file = glue('figs/{gene}_cc.pdf'), width=4,height=2)
}




