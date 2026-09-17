#--------------------------------
## prepare domain-specific deseq object
#--------------------------------
vsd <- varianceStabilizingTransformation(domain_dds, blind = F)
meanSdPlot(assay(vsd))

source('src/R/get_pc_plot.R')
p1 <- get_pc_plot(vsd, 'donor_id')
p2 <- get_pc_plot(vsd, 'tissue')
p3 <- get_pc_plot(vsd, 'age')

combined_plot <- wrap_plots(list(p3, p1, p2), ncol=3) + 
  plot_layout(guides = "keep")

ggsave(combined_plot, 
       filename = file.path('figs', glue('domain_counts_PCA.pdf')),
       width=12,height=4)