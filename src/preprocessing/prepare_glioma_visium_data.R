source('src/argo/00-load-packages.R')
library(Seurat)
library(UCell)
library(RColorBrewer)

# henrik <- read_rds('~/Desktop/brain-seurat-objects/heinrik.rds')
# idhm <- read_rds('~/Desktop/brain-seurat-objects/IDHm.rds')
# gbm <- read_rds('~/Desktop/brain-seurat-objects/GBM.rds')
# 
# henrik$subtype <- 'GBM'
# gbm$subtype <- 'GBM'
# idhm$subtype <- 'IDHmutant'
# 
# brain <- merge(henrik, c(gbm, idhm))
# brain <- JoinLayers(brain)


brain <- read_rds('~/Desktop/brain-seurat-objects/visium.rds')
#-------------------------------------------
#-------------------------------------------
# brain <- NormalizeData(brain)
# brain <- FindVariableFeatures(brain)
# brain <- ScaleData(brain)
# brain <- RunPCA(brain)
# brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
# brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
# DimPlot(brain, group.by = 'subtype')

#-------------------------------------------
## add ivy marker annotation
#-------------------------------------------
ivy_files <- grep('_markers.csv',list.files('data/sheets/'), value = T)
markers <- lapply(ivy_files, function(f) vroom(file.path('data/sheets', f), delim = ',')[[1]])
names(markers) <- c('CT', 'CTmvp', 'CTpan', 'LE')
markers <- lapply(markers, intersect, rownames(brain))
sapply(markers, length)

## score for marker genes with alternate methods, without counts imputation
brain <- AddModuleScore(brain, markers, name = 'ivy',
                              slot = 'data', ncores = 12)


FeaturePlot(brain, features = names(markers), ncol=2)


brain$ivy_niche <- as.factor(names(markers)[apply(brain@meta.data[paste0('ivy', 1:4)], 1, which.max)])
table(brain$ivy_niche, brain$ivy_region)
table(brain$ivy_niche, brain$subtype)

DimPlot(brain, group.by = 'ivy_niche')
SpatialDimPlot(brain, image.alpha = 1, images = 'slice1.5',
               group.by = 'ivy_niche', pt.size.factor = 2)
VlnPlot(brain, paste0('ivy', 1:4), pt.size = 0, group.by = 'subtype', ncol=2)
#-------------------------------------------
# don't need the scaled counts
#-------------------------------------------
brain@assays$Spatial@layers$scale.data <- NULL
write_rds(brain, '~/Desktop/brain-seurat-objects/visium.rds')
#-------------------------------------------
# bin data
#-------------------------------------------
source('src/R/bin_kdtree_adaptive.R')

ivy <- levels(brain$ivy_niche)
## clamp the enrichment scores at zero
xx <- pmax(as.matrix(brain@meta.data %>% dplyr::select(all_of(ivy))), 0)

## bin the data 
brain$bins <- bin_kdtree_adaptive(xx, target_spot_count = 5000)
(N <- length(levels(brain$bins)))

bin_cols <- colorRampPalette(brewer.pal(9, "Set3"))(N)
df <- xx %>% bind_cols(brain@meta.data %>% dplyr::select(
                                        ivy_niche, annotation1,
                                        condition = subtype, 
                                        sample_id = orig.ident,
                                        bin = bins))  %>% 
                         mutate(grid_key = paste(sample_id, bin, sep = '___'), 
                                condition_key = paste(condition, bin, sep = '___')) %>%
                         mutate(sample_id = factor(sample_id, 
                                                   levels = gtools::mixedsort(unique(sample_id))),
                                bin = factor(bin, levels = sort(unique(bin))), 
                                grid_key = factor(grid_key, 
                                                  levels = gtools::mixedsort(unique(grid_key))),
                                condition_key = factor(condition_key, 
                                                       levels = gtools::mixedsort(unique(condition_key)))) %>%
  arrange(grid_key)
                        


plotly::plot_ly(data = df, 
                x = ~CT, 
                y = ~CTpan, 
                z = ~LE, 
                color = ~condition_key, 
                colors = bin_cols,
                type = "scatter3d", 
                marker = list(size = 2.5, alpha=.5),
                mode = "markers")

# SpatialDimPlot(brain, group.by = 'bins', image.alpha = 0,
#                images = 'slice1', pt.size.factor = 3) + NoLegend()
# VlnPlot(brain, ivy, group.by = 'bins', pt.size = 0, ncol=2)
# VlnPlot(brain, 'nFeature_Spatial', group.by = 'condition', pt.size = 0)
# pheatmap(table(brain$bins, brain$ivy_niche))
# pheatmap(table(brain$bins, brain$subtype), name = 'number of spots')
# pheatmap(scale(table(brain$bins, brain$annotation1)))

#============================================================
## prepare environment for running argo
#============================================================

coords <- df %>% 
  group_by(bin) %>% 
  summarise(across(all_of(ivy), \(x) mean(x, na.rm = T)))

counts <- df %>%
  dplyr::count(grid_key, 
              MP = annotation1, name = 'n_cells') %>% 
  group_by(grid_key) %>%
  mutate(total_cells = sum(n_cells),
         proportion = n_cells / total_cells) %>%
  ungroup()

prop <- counts %>% select(grid_key, MP, proportion) %>%
  arrange(grid_key) %>% 
  pivot_wider(
    names_from = MP,           # Cell types → column names
    values_from = proportion,        # Proportions → cell values
    values_fill = 0,                 # NA → 0 (missing cell types)
    names_sort = TRUE                # Alphabetize cell type columns
  ) %>% 
  column_to_rownames('grid_key') %>%
  as.matrix()



bin_meta <- left_join(left_join(coords, df %>% 
  distinct(condition, sample_id, bin, grid_key, condition_key), 
  by = 'bin'), counts %>% distinct(grid_key, total_cells), 
  by = 'grid_key') %>% arrange(grid_key)

stopifnot('binned counts does not match with metadata' = all(rownames(prop) == bin_meta$grid_key))

## now filter out bins with very few cells
# hist(bin_meta$total_cells)
# quantile(bin_meta$total_cells, c(.1,.2, .25,.5))
# sum(bin_meta$total_cells < 5)

# remove binxslide combinations with fewer than 10 cells
bin_meta <- bin_meta %>% filter(total_cells > 10)
prop <- prop[bin_meta$grid_key, ]
#========================================================================
## build the full unobserved regular spatial grid (n_slides x N)
#========================================================================
full_grid <- expand_grid(
  sample_id =  levels(bin_meta$sample_id), ## first argument varies slowest
  bin =  levels(bin_meta$bin)) %>% 
  mutate(grid_key = paste(sample_id, bin, sep = '___'))

## make the same variable with bin_meta
bin_meta <- bin_meta %>% 
  arrange(match(grid_key, full_grid$grid_key))

## construct the permutation matrix projection_matrix
matched_indices <- match(bin_meta$grid_key, full_grid$grid_key)

n_obs <- nrow(bin_meta); n_full <- nrow(full_grid)
projection_matrix <- matrix(0, nrow = n_obs, ncol = n_full)
for(i in seq_len(n_obs)){
  projection_matrix[i, matched_indices[i]] <- 1
}
# image(projection_matrix)
#=========================================
## prepare the coregionalization matrix
#=========================================
slide_df <- bin_meta %>% 
  distinct(sample_id, condition)

W <- model.matrix(~0 + slide_df$condition)
B <- tcrossprod(W) 
image(B)

noise_diag <- as.numeric(as.factor(slide_df %>% 
                           pull(condition)))

## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + bin_meta$sample_id)
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)

stopifnot('counts matrix is not aligned with metadata' = all(rownames(prop) == bin_meta$grid_key))
#============================================================
# center the abundance ata
#============================================================
prop <- scale(prop, center = T, scale = F)

plot(prop[, 'cMESHyp']~ bin_meta$condition_key)
plot(prop[, 'cMESInf']~ bin_meta$condition_key)
plot(prop[, 'cMESVasc']~ bin_meta$condition_key)
plot(prop[, 'Neuron']~ bin_meta$condition_key)

#============================================================
mean_abx <- rowsum(prop, bin_meta$condition_key) %>%
  as.data.frame() %>% 
  rownames_to_column('condition_key') 
mean_abx <- pivot_longer(mean_abx, cols = -condition_key, 
                         names_to = 'cell_type', 
                         values_to = 'proportion') %>% 
  left_join(bin_meta %>% select(condition, condition_key, all_of(ivy)),
            by = 'condition_key', multiple = 'first')


pal43 <- colorRampPalette(brewer.pal(11, "Set3"))(43)
p <- ggplot(mean_abx, aes(x = CTmvp, 
                          y = proportion,  # Or n_cells
                          color = cell_type, 
                          group = cell_type)) +
  geom_line(linewidth = .5, alpha = 0.75, linetype=1) +
  # geom_point(size = 1.25, alpha = .75) +
  facet_wrap(~ interaction(condition),  # Or facet_wrap(~ age + batch)
             scales = "free_y", ncol=1) +
  # scale_y_continuous(labels = scales::percent) +  # If proportions
  labs(title = "Mean Cell Type Proportions vs. Isodepth Bins by Age",
       x = "Binned Isodepth", y = "Centered Proportion of Cells",
       color = "Cell Type") +
  scale_color_manual(values = pal43) +
  theme_classic() +
  guides(color = guide_legend(ncol = 12)) +  # Set number of columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p)

#============================================================
#the coordinate axis for input
isodepth_ <- bin_meta %>% select(all_of(ivy)) %>% distinct() %>% as.matrix()
#============================================================
# create environment to pass these variables for fitting
#============================================================
data_env <- new.env()
data_env$counts <- prop
data_env$metadata <- bin_meta
data_env$isodepth_ <- isodepth_
data_env$N <- N
data_env$projection <- projection_matrix
data_env$W <- W
data_env$B <- B
data_env$donor_identity <- cov_donor
data_env$noise_diag <- noise_diag
data_env$KERNEL = 'RBF'
data_env$FIT_DONOR <- FALSE
data_env$FIT_CONDITIONAL_NOISE <- TRUE
data_env$LAMBDA <- 1 # regularization for noise terms in GP fits

data_env$CONTRAST_VECTOR <- c(-1, 1)
data_env$CONTRAST <- 'mutant_vs_GBM'
data_env$RETURN_CONF_BANDS = TRUE # compute conf bands only if I am interested in recovering intervals
data_env$conf_level_alpha = 0.95
data_env$n_monte_carlo_simulations = 1000 # no.of monte carlo samples
data_env$features_to_fit <- colnames(prop)
data_env$test_coords <- NULL
data_env$PLOT_POSTERIOR <- FALSE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$RETURN_FITS <- TRUE
data_env$kernel_params <- NULL # if my parameters are already fit, I don't want to fit likelihood again

## save variables to memory
write_rds(data_env, file = 'output/visium_data_env.rds')
#============================================================
# data_env <- read_rds('output/visium_data_env.rds')
#============================================================
## source implementation modules
#============================================================
source('src/argo/fit_ARGO.R')
source('src/R/nll_and_gradients.R')
source('src/R/posterior_predictive_plots.R')
source('src/R/contrast-test-with-funGP.R')
source('src/R/generate_cc_matrix.R')
source('src/R/compute_analytical_pvalue.R')
source('src/R/derive_KLT_CI_bands.R')
source('src/R/compute_exceedence_statistics.R')
source('src/R/plot_conf_bands.R')

res.ARGO <- fit_ARGO(data_env)
write_rds(res.ARGO, file = 'output/mutant_vs_WT_results.rds')

# I don't need to fit the model within each slide separately.
# I can just jointly model across all slides and split by condition -- not modeling donor effects

