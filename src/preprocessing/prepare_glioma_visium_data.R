#' AddModuleScore for ivy niche enrichment scores
#' add a soft-plus transform, rather than clipping scores at zero, to retain 
#' identifiability across bins
#' scale the inputs before generating bins
#' get pooled abundance data by condition 
#' filter out bins with very few spots -- thresholding at 100 spots here.
#' define the projection and donor_covariance matrices.
#' scale abundance data before fitting. 
#' # I am running diff tests b/w conditions; pooling across all slides per condition
# have no use for condition_key here -- all information is summarized at level of condition_key
# coords per bin per condition; counts per bin per condition

source('src/argo/00-load-packages.R')
# source('src/R/bin_kdtree_adaptive.R')
set.seed(42)

brain <- read_rds('~/Desktop/brain-seurat-objects/visium.rds')
ivy <- levels(brain$ivy_niche)
## convert to relative scores for spots
raw_scores <- as.matrix(brain@meta.data[, ivy])
# Temperature tau = 0.15 controls sharpness toward the 4 vertices
beta <- 0.05
scores_softplus <- scale(beta * log1p(exp(raw_scores/beta)))
N_BINS <- 65
MIN_CELLS <- 30L # mask bins coming from with fewer than so many cells per slide

## bin with k-means voronoi cells 
message('now splitting data into bins:')
km <- kmeans(scores_softplus, centers = N_BINS, nstart = 25)
brain$bins <- factor(km$cluster)
#============================================================
df <- scores_softplus %>% 
  bind_cols(brain@meta.data %>% dplyr::select(
    ivy_niche, 
    MP = annotation1,
    condition = subtype, 
    sample_id = orig.ident,
  ))  %>% 
  mutate(condition = factor(condition),
         bin = factor(km$cluster)) %>%
  mutate(sample_id = factor(sample_id, 
                            levels = gtools::mixedsort(unique(sample_id))),
         bin = factor(bin, levels = sort(unique(bin)))) %>%
  rownames_to_column('spot_id') %>% 
  as_tibble()

# levels(df$condition) # "GBM"       "IDHmutant"
# levels(df$sample_id)
#============================================================
## prepare abundance data
#============================================================
message('now preparing abundance data:')
# define proportions relative to total number of cells in each bin; across all cell types. 
counts <- df %>%
  dplyr::count(condition, sample_id, bin, MP, name = 'n_cells') %>%
  # Use nesting() so complete() doesn't falsely assign slides to the wrong condition
  tidyr::complete(tidyr::nesting(condition, sample_id), bin, MP, fill = list(n_cells = 0)) %>%
  group_by(condition, sample_id, bin) %>%
  mutate(
    N_sb = sum(n_cells),
    supported = N_sb >= MIN_CELLS,
    
    # Within-slide, within-bin composition
    proportion_raw = n_cells / na_if(N_sb, 0),
    
    # Mask the whole slide–bin, not individual low-count MPs
    proportion = if_else(supported, proportion_raw, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    condition_key = paste(condition, bin, sep = "___"),
    grid_key = paste(sample_id, bin, sep = "___")
  ) %>%
  mutate(condition_key = factor(condition_key, 
                                levels = gtools::mixedsort(unique(condition_key))),
         grid_key = factor(grid_key, 
                           levels = gtools::mixedsort(unique(grid_key)))
  )

#============================================================
# sum(is.na(counts$proportion))
# plot(counts$proportion_raw~ counts$N_sb)
# plot(counts$proportion~ counts$N_sb)
# filter(counts, proportion == 1) %>% dplyr::count(bin, MP)
#============================================================
composition_check <- counts %>%
  filter(supported) %>%
  group_by(condition, sample_id, bin) %>%
  summarise(total = sum(proportion), .groups = "drop")

stopifnot(all(abs(composition_check$total - 1) < 1e-10))

# need to filter out spots with low counts per condition -- proportions are not valid; 
# or compute proportion by bin, not by conditino_key
prop <- counts %>% 
  select(grid_key, MP, proportion) %>%
  arrange(grid_key) %>% 
  pivot_wider(
    names_from = MP,           # Cell types → column names
    values_from = proportion,        # Proportions → cell values
    values_fill = 0,                 # NA → 0 (missing cell types)
    names_sort = TRUE                # Alphabetize cell type columns
  ) %>% 
  column_to_rownames('grid_key') %>%
  as.matrix()

#============================================================
# build aligned bin metadata
coords <- df %>% 
  group_by(bin) %>% 
  summarise(across(all_of(ivy), \(x) mean(x, na.rm = T)))

bin_meta <- counts %>% distinct(condition, sample_id, bin, grid_key, condition_key, N_sb) %>% 
  left_join(coords, by = 'bin') %>%
  arrange(grid_key)
#============================================================
#============================================================
# Verify means across all bins:
# tapply(prop[, 'cMESHyp'], bin_meta$condition, mean, na.rm=T)
# tapply(prop[, 'Vasc'],    bin_meta$condition, mean, na.rm=T)
# tapply(prop[, 'cOXPHOS'], bin_meta$condition, mean, na.rm=T)
#============================================================
prop_centered <- scale(prop, center = TRUE, scale = TRUE)
stopifnot('binned counts does not match with metadata' = all(rownames(prop_centered) == bin_meta$grid_key))
#============================================================
#============================================================
(N <- length(levels(bin_meta$bin)))
#============================================================
# drop bins with few counts
keep <- which(bin_meta$N_sb > MIN_CELLS)
prop_centered <- prop_centered[keep, ]
bin_meta <- bin_meta[keep, ]
stopifnot('check if all sample ids are in bin meta:' = all(levels(df$sample_id) %in% bin_meta$sample_id))
stopifnot('check if all bins are in bin meta:' = all(levels(df$bin) %in% bin_meta$bin))
#========================================================================
## build the full unobserved regular spatial grid (n_slides x N)
#========================================================================
message('preparing projection matrix')
full_grid <- expand_grid(
  sample_id =  levels(bin_meta$sample_id), ## first argument varies slowest
  bin =  levels(bin_meta$bin)) %>% 
  mutate(grid_key = paste(sample_id, bin, sep = '___')) %>% 
  mutate(grid_key = factor(grid_key, 
                           levels = gtools::mixedsort(unique(grid_key))))


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
dim(projection_matrix)
#=========================================
## prepare the coregionalization matrix
#=========================================
slide_df <- bin_meta %>% 
  distinct(sample_id, condition)

W <- model.matrix(~0 + slide_df$condition)
B <- tcrossprod(W) 
image(B)

noise_diag <- as.numeric(slide_df %>% 
                           pull(condition))

## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + bin_meta$sample_id)
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)
dim(cov_donor)

stopifnot('counts matrix is not aligned with metadata' = all(rownames(prop_centered) == bin_meta$grid_key))
#============================================================
#============================================================
#the coordinate axis for input
isodepth_ <- bin_meta %>% select(all_of(ivy)) %>% distinct() %>% as.matrix()
#============================================================
# create environment to pass these variables for fitting
#============================================================
data_env <- new.env()
data_env$counts <- prop_centered
data_env$metadata <- bin_meta
data_env$isodepth_ <- isodepth_
data_env$N <- N
data_env$projection <- projection_matrix
data_env$W <- W
data_env$B <- B
data_env$cov_donor <- cov_donor
data_env$noise_diag <- noise_diag # mapping from individual samples to condition
data_env$KERNEL = 'RBF'
data_env$FIT_DONOR <- TRUE
data_env$FIT_CONDITIONAL_NOISE <- TRUE
data_env$LAMBDA <- 1 # regularization for noise terms in GP fits

data_env$CONTRAST_VECTOR <- c(1, -1)
data_env$CONTRAST <- 'WT_vs_Mutant'
data_env$RETURN_CONF_BANDS = TRUE # compute conf bands only if I am interested in recovering intervals
data_env$conf_level_alpha = 0.95
data_env$n_monte_carlo_simulations = 3000 # no.of monte carlo samples
data_env$features_to_fit <- colnames(prop)
data_env$test_coords <- NULL
data_env$PLOT_POSTERIOR <- FALSE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$GET_POSTERIOR <- FALSE
data_env$RETURN_FITS <- TRUE
data_env$kernel_params <- NULL # if my parameters are already fit, I don't want to fit likelihood again

## save variables to memory
write_rds(data_env, file = glue('output/glioma_kmeans_bysample_N_{N_BINS}_data_env.rds'))
#============================================================

#============================================================
## source implementation modules
#============================================================
res.ARGO <- fit_ARGO(data_env)
print(res.ARGO$params)
write_rds(res.ARGO, glue('output/glioma_kmeans_bysample_N_{N_BINS}_results.rds'))
#============================================================
# data_env <- read_rds('output/glioma_condition_data_env.rds')
# list2env(as.list(data_env), envir = globalenv())
# res.ARGO <- read_rds('output/WT_vs_mutant_condition_results.rds')
#=================================================
## visualize results
#=================================================
funGP_res <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'feature') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation)) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 
funGP_res
#=================================================
## where the significant deviation lies
#=================================================
p <- ggplot(funGP_res %>% arrange(desc(stats)), 
            aes(x = Rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))},
             color = "red", size = 2.5, alpha=0.5) +
  
  # Add the non-overlapping labels
  geom_text_repel(
    data = . %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))},
    aes(label = feature),
    size = 5, 
    family = 'Arial',
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "red", 
    max.overlaps = 50
  ) +
  # Labels and clean styling
  labs(
    x = "Rank",
    y = "Signed neglogp",
    title = glue("ARGO Rank Statistics")
  ) +
  theme_classic(base_family = 'Arial') +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size=14)
  )

print(p)
