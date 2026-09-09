pvalsDf <- pvalsDf %>% 
  mutate(baseline_abx = map2_dbl(
    pvalsDf$index, 
    pvalsDf$cell, 
    ~ baseline_abx[as.integer(.x), as.character(.y)]
  )) 


cx <- 'adult'
roi <- idx_upreg$`MSN-D1-1`
dx <- dxlist[[2]]
cx1 <- Cells(dx)[dx$clust_annot == cell]
cx2 <- Cells(dx)[dx$binned_isodepth %in% roi]

p1 <- DimPlot(dx, reduction = 'X_spatial_coords', 
        group.by = 'clust_annot', cells.highlight = cx1) +
  ggtitle('') +
  labs(x = "", y = "", color = cell)

p2 <- DimPlot(dx, reduction = 'X_spatial_coords', 
        group.by = 'clust_annot', cells.highlight = cx2) + 
  ggtitle('') +
  labs(x = "", y = "", color = "ROI")

pdf(glue('figs/{cx}_{cell}_dist.pdf') ,width = 12, height = 6)
wrap_plots(p1+p2) + plot_annotation(title = glue('{cx} brain: {cell}'))
dev.off()
