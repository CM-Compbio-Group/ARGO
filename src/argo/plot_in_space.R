#===================================================================
## visualize some representative genes in space
#===================================================================
genes_to_plot <- c('Serpina3n', 'Irf7', 'Foxp2', 'Sox10')
library(Seurat)
dxlist <- read_rds('~/Desktop/slide_list.rds')
dxlist <- dxlist[c(3,4,5,6)]
dxlist <- dxlist[c(1,3)]

plot_domain <- function(dx){
  p <- DimPlot(dx,  group.by  = 'tissue',
               reduction = 'X_spatial_coords', pt.size=.5 
  ) + 
    ggtitle('') + 
    labs(x = '', y = '') + 
    theme_minimal(base_family = 'Arial') + 
    theme(axis.text = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(p)
}

plots <- lapply(dxlist, plot_domain)
row_plot <- wrap_plots(plots, ncol = 2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = 'Tissue Annotation')

pdf(glue('figs/spatial_domains.pdf'), 
    width = 8, height = 3.5)
plot(row_plot)
dev.off()

tc <- do.call('rbind', lapply(dxlist, function(dx) Embeddings(dx, reduction = 'X_spatial_coords')))
colnames(tc) <- c('x_', 'y_')

dfm <- do.call(rbind, lapply(dxlist, function(dx) dx@meta.data))
dfm <- dfm %>% 
  mutate(coord_x = filter(metadata, cell_id %in% rownames(dfm)) %>%
           dplyr::slice(match(rownames(dfm), cell_id)) %>% 
           pull(coord_x)
  )
           

alter_seurat_meta <- function(dx){
  dx@meta.data <- dx@meta.data %>% 
    mutate(coord_x = metadata %>% filter(cell_id %in% Cells(dx)) %>% 
             dplyr::slice(match(rownames(dfm), cell_id)) %>% pull(coord_x))
  dx
}

dxlist <- lapply(dxlist, alter_seurat_meta)


FeaturePlot(dxlist[[2]], feature = 'coord_x',
            reduction = 'X_spatial_coords', 
            alpha = .5,
            pt.size=.5) + 
  scale_color_gradientn(
    colors = c("lightgray", 'red1', "red2",'red3', "darkred"),
  )  

## prepare common metadata file
dfm <- do.call(rbind, lapply(dxlist, function(dx) dx@meta.data))

## get roi for each gene
roi <- lapply(funGP_fits[genes_to_plot], function(kxx) isodepth_[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)

## extract counts for genes of interest across all slides
cts <- do.call('cbind', lapply(dxlist, function(dx) LayerData(dx)[genes_to_plot, ]))
coord_x <- dfm$coord_x
# cdx <- paste(dfm$donor_id, as.numeric(as.factor(dfm$coord_x)), sep = '_')

## need to set thresholds across all three age groups by roi
highlight_roi <- function(dx, rx, gene_thresh, gene){
  expr <- LayerData(dx)[gene, ]
  alpha_map <- ifelse(dx$coord_x %in% rx & expr, .75, 5e-3)
  size_map <-  ifelse(dx$coord_x %in% rx, .5, 1e-2)
  p <- FeaturePlot(dx, feature = gene,
                   reduction = 'X_spatial_coords', 
                   # min.cutoff = 0,
                   # max.cutoff = 2,
                   alpha = alpha_map,
                   pt.size=size_map) + 
    scale_color_gradientn(
      colors = c("lightgray",  "red2", "darkred"),
      limits = gene_thresh[c(1,3)],
      breaks = gene_thresh
    ) + 
    ggtitle('') + 
    labs(x = '', y = '') + 
    theme_minimal(base_family = 'Arial') + 
    theme(axis.text = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}

figs <- lapply(genes_to_plot, function(gene) {
  message(gene)
  rx <- roi[[gene]];
  cx <- cts[gene, (coord_x %in% rx)];
  gene_thresh <- quantile(cx[cx > 0], c(0.5, .75, 0.99));
  lapply(dxlist, highlight_roi, rx, gene_thresh, gene )
}
)

names(figs) <- genes_to_plot
plots <- lapply(genes_to_plot, function(gene) {
  message(gene)
  row_plot <- wrap_plots(figs[[gene]], ncol = 4) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{gene}'))
  wrap_elements(row_plot)
}
)

dynamic_height <- max(6, length(genes_to_plot) * 3.5)
pdf(glue('figs/{CONTRAST}_partial_extent_spatial_ROI_quantiles.pdf'), 
    width = 14, height = dynamic_height)
print(wrap_plots(plots, nrow = length(genes_to_plot)) + 
        plot_layout(guides = "collect") + 
        plot_annotation(title = glue('{CONTRAST}')))
dev.off()


# boxplot(cts['Ctps1', ]~ cdx)

dx <- dxlist[[5]]
expr <- LayerData(dx)[gene, ]
alpha_map <- ifelse(expr, .95, 1e-1)
FeaturePlot(dx, feature = gene,
            reduction = 'X_spatial_coords', 
            pt.size=.5, alpha = alpha_map,
) +
  scale_color_gradientn(
    colors = c("lightgray",  "red2", "darkred"),
  ) + 
  ggtitle('old mouse') + 
  labs(x = '', y = '') + 
  theme_minimal(base_family = 'Arial') + 
  theme(axis.text = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ad <- merge(dxlist[[1]], dxlist[2:6])
ad <- JoinLayers(ad)
VlnPlot(ad, 'Ctps1', group.by = 'donor_id', pt.size=0)
DotPlot(ad, genes_to_plot, cols = c('lightgrey', 'darkred'), group.by = 'donor_id')

tapply(LayerData(ad)['Ctps1', ], ad$donor_id, function(x){
  c(
    pct_nonzero = round(mean(x > 0) * 100, 1),
    mean_overall = round(mean(x), 3),
    mean_nonzero = round(mean(x[x > 0]), 3)
  )
})

table(ad$donor_id, ad$age)



library(ggplot2)
library(viridis)

# Bind Ctps1 expression to your spatial metadata
plot_df <- dfm %>%
  bind_cols(t(as.matrix(cts)), tc)

# Option A: Color by Local Mean Expression (what pseudobulk/ARGO sees)
ggplot(plot_df, aes(x = x_, y = y_, z = Foxp2)) +
  stat_summary_hex(fun = mean, bins = 75) +
  scale_fill_viridis_c(option = "viridis", name = "Mean\n Foxp2", 
                       limits = c(0, 2)) +
  facet_wrap(~ age + donor_id) +
  coord_equal() +
  theme_void(base_size = 14)

# Option B: Color by Local % Expressing Cells (if the shift is driven by dropout rate)
ggplot(plot_df, aes(x = x_, y = y_, z = as.numeric(Serpina3n > 0))) +
  stat_summary_hex(fun = mean, bins = 60) +
  scale_fill_viridis_c(option = "inferno", labels = scales::percent, name = "% Serpina3n+") +
  facet_wrap(~ age + donor_id) +
  coord_equal() +
  theme_void(base_size = 14)


ggplot(plot_df, aes(x = x_, y = y_, color=coord_x)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  facet_wrap(~ age + donor_id) +
  coord_equal() +
  theme_void()
