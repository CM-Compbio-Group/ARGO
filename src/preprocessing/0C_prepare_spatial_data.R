## prepare spatial data
library(anndataR)
library(Seurat)

## read raw h5ad file and convert to Seurat object
# ad <- read_h5ad('data/raw/integrated_anndata.h5ad')
ad <- read_h5ad('~/Desktop/integrated_anndata.h5ad')
ad <- ad$as_Seurat()

DimPlot(ad, reduction = 'X_umap', group.by = 'donor_id')
DimPlot(ad, reduction = 'X_spatial_coords', group.by = 'donor_id')

levels(ad$batch) <- paste('slide', levels(ad$batch), sep = '-')
table(ad$batch, ad$donor_id)

## visualize domains in each slide
figs <- list()
for(bx in levels(ad$batch)){
  figs[[bx]] <- DimPlot(subset(ad, batch == bx), reduction = 'X_spatial_coords', group.by = 'tissue') + 
    ggtitle(glue('slide_{bx}')) + 
    labs(x = '', y = '') + 
    theme_minimal(base_family = 'Arial') + 
    theme(axis.text = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

image_plots <- wrap_plots(figs, ncol = 5) + 
  plot_layout(guides = "collect") 

pdf(glue('figs/spatial_domains_by_slide.pdf'), width = 15, height = 15)
print(image_plots)
dev.off()

#-------------------------------------
## extract slides of interest
#-------------------------------------
table(ad$unique_slice, ad$age)
# visualize domain compositions of all slides
plots <- list()
for(id in unique(ad$unique_slice)){
  message(id)
  dx <- subset(ad, subset = unique_slice == id)
  plots[[id]] <- DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'tissue') 
}

pdf('figs/slice_tissues.pdf', width = 15, height = 15)
wrap_plots(plots, ncol=5, guides = "collect")
dev.off()

## identify tissues to use for visualization
idxlist <- list('Mouse-4_2', 'Mouse-7_2', 'Mouse-10_0', 'Mouse-12_0', 'Mouse-5_0', 'Mouse-6_1')

dxlist <- lapply(idxlist, function(id) subset(ad, subset = unique_slice == id))

# dxlist <- lapply(dxlist, function(dx) {Idents(dx) <- 'cell_type'; dx})
plots <- lapply(dxlist, function(dx) DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'tissue') + 
                  ggtitle(glue('donor_{dx$donor_id[1]}_{dx$batch[1]}')) + 
                  labs(x = '', y = '') + 
                  theme_minimal(base_family = 'Arial') + 
                  theme(axis.text = element_blank(), 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank()))

p <- wrap_plots(plots, ncol = 2) + plot_layout(guides = "collect")

pdf('figs/chosen_slides.pdf', width = 8, height = 10)
print(p)
dev.off()

## save this list of slides to file
write_rds(dxlist, 'data/processed/slide_list.rds')
write_rds(ad@meta.data, 'data/processed/ad_metadata.rds')
#-----------------------------------------------------------
## learn the delta parameter to use for mapping test points to spots
#-----------------------------------------------------------
isx <- sort(unique(ad$binned_isodepth))
(delta <- quantile(diff(isx), .25))

#-------------------------------------
## plot cell types of interest
#-------------------------------------
dx <- dxlist[[1]]
## see how the x_test isodepth values compare with binned isodepths
FeaturePlot(dx, reduction = 'X_spatial_coords', features = 'binned_isodepth')
hist(dx$binned_isodepth)

## plot for specific celltype
roi <- idx_downreg$OPC

## check delta b/w bins
isx <- sort(unique(dx$binned_isodepth))
## set delta to twice dist b/w adjacent bins
(delta <- 2*min(diff(isx)))

## highlight spots w/in delta of roi
## only need to find which bins are within delta of roi
disx.roi <- abs(outer(isx, roi, FUN = '-')) < delta
uniq.isx <- isx[which(apply(disx.roi, 1, any))]

## highlight only these cells
dx$inROI <- ifelse(dx$binned_isodepth %in% uniq.isx, T, F)
DimPlot(dx, reduction = 'X_spatial_coords', group.by = c('tissue', 'inROI'))

plots <- lapply(dxlist, plot_cell_states, 'oligodendrocyte')
wrap_plots(plots, ncol=3, guides = 'collect')
