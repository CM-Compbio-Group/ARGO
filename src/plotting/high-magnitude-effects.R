## filter to genes of interest
df <- filter(pvalsDf, net_poe < .5 & abs(deviation) > 1 &  neglogp > -log(.05))

pdf(glue('figs/{CONTRAST}_FCmap_select.pdf'), width = 6, height = 4)
Heatmap(stat_lfc[df$gene, ], cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=4))
dev.off()

#-----------------------------------------------
genes_to_plot <- df$gene
condition <- 'discordant_trends'
length(genes_to_plot)

source('src/scripts/02_visualize_fitted_GPs.R')
px <- wrap_plots(plots, ncol = 5) + 
  plot_layout(guides = "collect")

ggsave(px,filename = file.path('figs', glue('{CONTRAST}_posterior_fit_{condition}.pdf')),
       width=18,height=8)

#-----------------------------------------------
conf_band_plots <- lapply(genes_to_plot, function(gene) 
  plot_CI_bands(gene, x = x_test, 
                band = exact_pvals[[gene]][['band']],
                stat = exact_pvals[[gene]][['statistic']], 
                alpha = 95)
)

py <- wrap_plots(conf_band_plots, ncol = 5) + 
  plot_layout(guides = "collect")

ggsave(py,filename = file.path('figs', glue('{CONTRAST}_CI_bands_{condition}.pdf')),
       width=18,height=8)

#-----------------------------------------------
#-----------------------------------------------
## get roi for each gene
roi <- lapply(exact_pvals[genes_to_plot], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)

## find what domains rois are mapping to
dfm <- rbind(dxlist[[1]]@meta.data, dxlist[[2]]@meta.data)
tt <- do.call('rbind', lapply(roi, function(rx) {
      filter(dfm, binned_isodepth %in% rx) %>%
        with(table(tissue))
    })
  )

tt <- tt / rowSums(tt)
pdf(glue('figs/{CONTRAST}_domain_distribution_select.pdf'), width = 6, height = 5)
Heatmap(tt, name = 'percent \nROI', 
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),)
dev.off()

#-----------------------------------------------
# assign genes to single-max domain
xt <- colnames(tt)[apply(tt, 1, function(t) order(t, decreasing=T)[2])]
names(xt) <- rownames(tt)

#-----------------------------------------------
## plot domain ROIs
#-----------------------------------------------
highlight_domain <- function(dx, rx, gene_thresh, gene){
  expr <- LayerData(dx)[gene, ]
  alpha_map <- ifelse(dx$tissue %in% rx & expr, .25, 1e-2)
  size_map <-  ifelse(dx$tissue %in% rx, .25, 1e-2)
  p <- FeaturePlot(dx, feature = gene,
                   reduction = 'X_spatial_coords', 
                   # min.cutoff = gene_thresh[1],
                   # max.cutoff = gene_thresh[2],
                   alpha = alpha_map,
                   pt.size=size_map) + 
    scale_color_gradientn(
      colors = c("lightgray", "red2", "darkred"),
      limits = gene_thresh[c(1,5)],
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
  rx <- xt[gene];
  cx <- cts[gene, domainx == rx];
  # plot(cx~ filter(dfm, tissue == rx) %>% pull(unique_slice) %>% droplevels())
  gene_thresh <- quantile(cx[cx > 0], c(.25, .5, .75, 0.90, .99));
  lapply(dxlist[3:6], highlight_domain, rx, gene_thresh, gene )
  }
)


names(figs) <- genes_to_plot
plots <- lapply(genes_to_plot, function(gene) {
  message(gene)
  row_plot <- wrap_plots(figs[[gene]], ncol = 4) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{gene} ROI'))
  # Protect the row assembly from being cleared by the top-level layout
  wrap_elements(row_plot)
}
)

pdf(glue('figs/{condition}_spatial_domain.pdf'), width = 12, height = 30)
print(wrap_plots(plots, nrow=10) + 
        plot_annotation(title = glue('{CONTRAST}')))
dev.off()


