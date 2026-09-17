## focus on mixed-signature genes
source('src/R/plot_conf_bands.R')

gx <- filter(df1, poe_up > 0 & poe_down > 0) %>% pull(gene)

pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_heatmap_mixed_signature.pdf'), width = 6, height = 4)
myc <- colorRamp2(breaks = c(-1, 0, 1), colors = c("blue3", "lightgray", "red3"))
Heatmap(stat_lfc[gx, ], 
        show_column_names = F,
        cluster_columns = F, name = 'FC \nstatistic',
        col = myc,
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=4))
dev.off()



conf_band_plots <- lapply(gx, function(gene) 
  plot_CI_bands(gene, x = x_test, 
                band = exact_pvals[[gene]][['band']],
                stat = exact_pvals[[gene]][['statistic']], 
                alpha = 95)
)

ggsave(wrap_plots(conf_band_plots, ncol = 3) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_conf_bands_mixed_sig.pdf')),
        width=12,height=8)

#------------------------------------------------------
#------------------------------------------------------
df1 %>% filter(neglogp > -log(.05) & abs(net_auec) > .25)