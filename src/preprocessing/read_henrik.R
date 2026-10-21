library(Seurat)
library(tidyverse)
library(vroom)
library(glue)

meta<- vroom('data/sheets/metaprogram_annotations_df.csv')
henrik_meta <- vroom('data/sheets/sampleID_conversion_GBM_Ravi22.csv')
id <- paste0('#UKF',str_sub(henrik_meta$SampleID_Ravi2022, start = 2, end=4), '_T')
henrik_meta$id <- id

#-------------------------------------------
#-------------------------------------------
data_dir <- '/Volumes/umms-wna/Collaborator Data_Freiburg/10XVisium/'

henrik.drx <- list.dirs(data_dir, full.names = F, recursive = F)
drx <- grep(pattern = paste(id, collapse = '|'), henrik.drx, value = T)
new_cell_id <- henrik_meta$SampleID_this_study

# read samples
slide.list <- list()
for(i in seq_along(drx)){
  message(new_cell_id[i])
  obj <- Load10X_Spatial(data.dir = file.path(data_dir, drx[i], 'outs'))
  obj$orig.ident <- new_cell_id[i]
  slide.list[[drx[i]]] <- obj
}

#-------------------------------------------
## map cells to their new identity
#-------------------------------------------
sapply(slide.list, ncol)

## merge into single object
brain <- merge(slide.list[[1]], slide.list[2:length(drx)], add.cell.ids = new_cell_id)
brain <- JoinLayers(brain)

VlnPlot(brain, 'nFeature_Spatial', group.by = 'orig.ident', pt.size = 0)
## append cell type annotations from alissa
cells <- intersect(Cells(brain), meta$full_id)
brain <- subset(brain, cells= cells)
all(Cells(brain) == cells)

brain@meta.data[c('annotation1', 'annotation2')] <- meta %>% 
                           dplyr::slice(match(cells, full_id)) %>%
                           select(annotation1, annotation2)

SpatialDimPlot(brain, group.by = 'annotation1', images = 'slice1.13',
               pt.size.factor = 2.5, image.alpha = 0)

write_rds(brain, 'data/temp/heinrik.rds')
