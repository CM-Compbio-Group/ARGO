## take diff testing results and extract rois of interest for each gene
## rois defined as regions falling outside the desired conf bands
#-------------------------------------
#-------------------------------------
#-------------------------------------
## read list of seurat objects with representative images
# dxlist <- read_rds('data/processed/slide_list.rds')
#-------------------------------------

## read diff testing results
pvalsDf <- pvalsDf %>% arrange(desc(stats), desc(effect_size))

upreg_cells <-  pvalsDf %>% pull(cell) %>% head(n = 5)
downreg_cells <- pvalsDf %>% pull(cell) %>% tail(n = 5)

## for these cells, identify rois to plot
idx_upreg <- lapply(exact_pvals[upreg_cells], function(kxx) which(kxx$statistic > kxx$band))
idx_downreg <- lapply(exact_pvals[downreg_cells], function(kxx) which(kxx$statistic < -kxx$band))

## get corresponding isodepth values
idx_upreg <- lapply(idx_upreg, function(ix) x_test[ix])
idx_downreg <- lapply(idx_downreg, function(ix) x_test[ix])

## how large are the rois
sapply(idx_upreg, length)
sapply(idx_downreg, length)

## plot these isodepth bins in slides
## for each celltype, get list of plots
upreg_figs <- lapply(names(idx_upreg), function(cell) 
{lapply(dxlist, plot_roi, idx_upreg[[cell]], cell)})

names(upreg_figs) <- names(idx_upreg)
figs <- lapply(names(upreg_figs), function(cell) {
  row_plot <- wrap_plots(upreg_figs[[cell]], ncol = 3) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{cell} ROIs'))
  
  # Protect the row assembly from being cleared by the top-level layout
  wrap_elements(row_plot)
}
)

pdf(glue('figs/{contrast}_top5_upreg.pdf'), width = 12, height = 20)
p <- wrap_plots(figs, nrow = 5) + 
  plot_annotation(title = glue('{contrast} upregulated'))
print(p)
dev.off()


## show the quantified effect sizes for each cell highlighted
pvalsDf %>% head(n = 5)
#-------------------------------------
## repeat for downreg celltypes
#-------------------------------------
downreg_figs <- lapply(names(idx_downreg), function(cell) 
{lapply(dxlist, plot_roi, idx_downreg[[cell]], cell)})

names(downreg_figs) <- names(idx_downreg)
figs <- lapply(names(downreg_figs), function(cell) {
  row_plot <- wrap_plots(downreg_figs[[cell]], ncol = 3) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{cell} ROIs'))
  
  # Protect the row assembly from being cleared by the top-level layout
  wrap_elements(row_plot)
}
)

pdf(glue('figs/{contrast}_top5_downreg.pdf'), width = 12, height = 20)
p <- wrap_plots(figs, nrow = 5) + 
  plot_annotation(title = glue('{contrast} downregulated'))
print(p)
dev.off()

## get plot for specific celltype
# plots <- lapply(dxlist, plot_roi, idx, cell = 'Olig-2')
# p <- wrap_plots(plots, ncol = 3) + 
#   plot_layout(guides = "collect") + 
#   plot_annotation(title = glue('{cell} ROIs'))
# pdf(glue('figs/{cell}_roi.pdf'), width=12, height=4)
# print(p)
# dev.off()

## get stats
pvalsDf %>% tail(n = 5) %>% select(cell, neglogp, deviation, baseline_abx, effect_size)


pvalsDf %>% filter(baseline_abx > .1)
