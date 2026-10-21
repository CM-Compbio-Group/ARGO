library(Seurat)
library(tidyverse)
library(vroom)
library(glue)

meta<- vroom('data/sheets/metaprogram_annotations_df.csv')

#-------------------------------------------
#-------------------------------------------
data_dir <- '/Volumes/umms-congma1/users/vravik/GBM-analysis/tirosh_cell_visium/data/raw/data_zenodo_RH/IDHmut/'

IDHm.drx <- list.dirs(data_dir, full.names = F, recursive = F)
# remove the LQ sample
IDHm.drx <- IDHm.drx[-24]

# read samples
slide.list <- list()
for(drx in IDHm.drx){
  message(drx)
  obj <- Load10X_Spatial(data.dir = file.path(data_dir, drx, 'outs'))
  obj$orig.ident <- drx
  slide.list[[drx]] <- obj
}

#-------------------------------------------
## map cells to their new identity
#-------------------------------------------
sapply(slide.list, ncol)

## merge into single object
brain <- merge(slide.list[[1]], slide.list[2:length(IDHm.drx)], add.cell.ids = IDHm.drx)
brain <- JoinLayers(brain)

VlnPlot(brain, 'nFeature_Spatial', group.by = 'orig.ident', pt.size = 0)

## append cell type annotations from alissa
cells <- intersect(Cells(brain), meta$full_id)
brain <- subset(brain, cells= cells)

brain@meta.data[c('annotation1', 'annotation2')] <- meta %>% 
  dplyr::slice(match(cells, full_id)) %>%
  select(annotation1, annotation2)

SpatialDimPlot(brain, group.by = 'annotation1', images = 'slice1',
               pt.size.factor = 2.5, image.alpha = 0)

SpatialDimPlot(brain, group.by = 'annotation1', images = 'slice1.2',
               pt.size.factor = 2.5, image.alpha = 0)

write_rds(brain, 'data/temp/IDHm.rds')
