filter(df1, gene %in% cas_positive)
df2 <- df2 %>% column_to_rownames('gene')
df2 <- df2[rownames(df1), ]
names(df1) <- paste('fungp', names(df1), sep='_')
names(df2) <- paste('deseq', names(df2), sep='_')


dfx <- cbind(df1, df2) %>% rownames_to_column('gene') %>% as_tibble()
plot(dfx$fungp_stats, dfx$deseq_stats)
plot(dfx$fungp_deviation,dfx$deseq_maxlogFC)

idx_dn <- which((dfx$fungp_stats < log(.05)) & (dfx$deseq_stats > log(.05)))
idx_dn_dev <- which(((dfx$fungp_stats < log(.05)) & (dfx$fungp_deviation < -1)) & (dfx$deseq_stats > log(.05)))

plot(dfx$deseq_stats[idx_dn], pch = 20, cex=1.25, 
     ylab = 'neglog combined pval')
genes[idx_dn]

plot(dfx$deseq_stats[idx_dn_dev], pch = 20, cex=1.25,
     ylab = 'neglog combined pval')

write.table(dfx[idx_dn, ] %>% pull(gene), 'output/files/fungp_dn_delta.txt', quote = F, row.names = F)

## plot the genes with opposing signals
idx_dn_up <- idx_dn[which(dfx[idx_dn, 'deseq_stats'] > -log(.05))]
View(dfx[idx_dn_up, ])

all(rownames(logFC) == dfx$gene)
pdf('figs/DEseq_deviant_LFC_heatmap.pdf', width = 6, height = 4)
Heatmap(logFC[idx_dn_up, ], name = 'LFC', 
        show_row_names = T,
        cluster_columns = F, cluster_rows = F,
        column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
        heatmap_legend_param = list(
          title = "LogFC",
          title_gp = gpar(fontsize = 16, fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
          legend_height = unit(4, "cm"), # Adjust height of the color bar
          grid_width = unit(5, "mm")      # Adjust width of the color bar
        ))
dev.off()


conf_band_plots <- lapply(genes[idx_dn_up], function(cell) 
  plot_CI_bands(cell, x = x_test, 
                band = exact_pvals[[cell]][['band']],
                stat = exact_pvals[[cell]][['statistic']], 
                alpha = 95)
)


ggsave( wrap_plots(conf_band_plots, ncol = 2) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{contrast}_conf_bands_deviant.pdf')),
        width=8,height=6)

## visualize fitted gps on raw data
plots <- list()
for(gene in genes[idx_dn_up]){
  message(gene)
  y <- bin_meta %>% pull(gene)
  ## pass corresponding gp parameters
  plots[[gene]]  <- posterior_predictive_plots(cell = gene, y = y,
                                               par = params[gene, ], 
                                               x_observed = x_observed,
                                               x_test = x_test,
                                               P_obs = P_obs,
                                               dsqx = dsqx,
                                               dsqt = dsqt,
                                               B = B, 
                                               N = N, 
                                               M = M,
                                               slides_per_group = slides_per_age,
                                               noise_diag = noise_diag,
                                               bin_meta = bin_meta)
}



ggsave( wrap_plots(plots, ncol = 2) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('gp_posterior_fits_deviant.pdf')),
        width=12,height=6)

ggsave(file.path('figs', glue('{gene}_gp_posterior_fit.pdf')),plot = p,
       width=8,height=4)
