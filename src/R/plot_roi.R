## passed a seurat object, highlight the roi and the cell type of interest
plot_roi <- function(dx, roi, cell){
  
  ## set color palette
  tx <- ifelse(dx$binned_isodepth %in% roi, T, F)
  palx <- rep('F3', ncol(dx))
  palx[(tx) & (dx$clust_annot == cell)] <- 'F1' # target cells in ROI
  palx[(tx) & (dx$clust_annot != cell)] <- 'F2' # non-target cells in ROI
  palx[!(tx)] <- 'F3' ## not in ROI
  dx$palx <- factor(palx)
  
  # palette <- RColorBrewer::brewer.pal(n = 4, 'Set2')
  palette <- c('red3', 'skyblue3', 'gray')
  names(palette) <- levels(palx)
  
  alpha_map <- recode(palx, "F1" = .75,
                      "F2" = .25,
                      "F3" = .1)
  
  p <- DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'palx',
               cols = palette, 
               alpha = alpha_map,
               # alpha = ifelse(palx %in% c('F1'), 0.75, 0.25),
               pt.size = alpha_map
               # pt.size = ifelse(palx %in% c('F1'), 0.75, 0.25)
               ) + 
    ggtitle('') + 
    labs(x = '', y = '') + 
    scale_color_manual(
      name = cell,
      values = palette, 
      labels = c('F1' = "target cell; in roi", 
                 'F2' = "non-target cell; in roi", 
                 'F3' = "outside roi"),
      guide = guide_legend(override.aes = list(alpha = 1, size = 3)) # Keeps dots legible
      
    )
  return(p)
}



