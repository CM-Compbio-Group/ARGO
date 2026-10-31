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
library(Seurat)
library(RColorBrewer)
set.seed(42)

brain <- read_rds('~/Desktop/brain-seurat-objects/visium.rds')
ivy <- levels(brain$ivy_niche)
## convert to relative scores for spots
raw_scores <- as.matrix(brain@meta.data[, ivy])
# Temperature tau = 0.15 controls sharpness toward the 4 vertices
beta <- 0.05
scores_softplus <- scale(beta * log1p(exp(raw_scores/beta)))
## bin with k-means voronoi cells 
km <- kmeans(scores_softplus, centers = 100, nstart = 25)
brain$bins <- factor(km$cluster)

#============================================================
df <- scores_softplus %>% 
  bind_cols(brain@meta.data %>% dplyr::select(
                                        ivy_niche, 
                                        MP = annotation1,
                                        condition = subtype, 
                                        sample_id = orig.ident,
                                        bin = bins))  %>% 
                         mutate(condition = factor(condition)) %>%
                         mutate(sample_id = factor(sample_id, 
                                            levels = gtools::mixedsort(unique(sample_id))),
                                bin = factor(bin, levels = sort(unique(bin)))) %>%
                        rownames_to_column('spot_id') %>% 
                        as_tibble()

levels(df$condition) # "GBM"       "IDHmutant"
#============================================================
## prepare abundance data
#============================================================
## account for relative abundance of WT and mutant spots
cond_weights <- df %>%
  dplyr::count(condition, name = "total_cond_spots") %>%
  mutate(cond_weight = total_cond_spots / sum(total_cond_spots))

# Total spots in each bin across BOTH conditions combined
bin_totals <- df %>%
  dplyr::count(bin, name = "total_bin_spots")

View(df %>% dplyr::count(condition, MP, name = 'n_cells'))

# define proportions relative to total number of cells in each bin; across all cell types. 
# number of cells from each condition in this bin
counts <- df %>%
  dplyr::count(condition, bin, MP, name = 'n_cells') %>%
  tidyr::complete(condition, bin, MP, fill = list(n_cells = 0)) %>%
  # 1. Local spot count of condition c in bin b (N_{c,b})
  group_by(condition, bin) %>%
  mutate(
    condition_key = paste(condition, bin, sep = '___'),
    N_cb = sum(n_cells)
  ) %>%
  ungroup() %>%
  # 2. Global spot count of condition c (N_c)
  group_by(condition) %>%
  mutate(
    N_c = sum(n_cells),
    state_purity = n_cells / pmax(N_cb, 1),
    cond_bin_density = N_cb / N_c
  ) %>%
  ungroup() %>%
  # 3. Total local density in bin b across both conditions: (N_{GBM,b}/N_GBM + N_{IDH,b}/N_IDH)
  group_by(bin) %>%
  mutate(
    # Divide by 17 (n_distinct(MP)) because cond_bin_density is repeated across the MP rows
    total_local_bin_density = sum(cond_bin_density) / n_distinct(MP),
    local_niche_share = cond_bin_density / total_local_bin_density,
    # Final local proportion in [0, 1] that sums to 1 inside every bin b
    proportion = state_purity * local_niche_share
  ) %>%
  ungroup() %>%
  mutate(
    condition_key = factor(condition_key, levels = gtools::mixedsort(unique(condition_key)))
  )

# need to filter out spots with low counts per condition -- proportions are not valid; 
# or compute proportion by bin, not by conditino_key
prop <- counts %>% 
  select(condition_key, MP, proportion) %>%
  arrange(condition_key) %>% 
  pivot_wider(
    names_from = MP,           # Cell types → column names
    values_from = proportion,        # Proportions → cell values
    values_fill = 0,                 # NA → 0 (missing cell types)
    names_sort = TRUE                # Alphabetize cell type columns
  ) %>% 
  column_to_rownames('condition_key') %>%
  as.matrix()

#============================================================
# build aligned bin metadata
coords <- df %>% 
  group_by(bin) %>% 
  summarise(across(all_of(ivy), \(x) mean(x, na.rm = T)))

bin_meta <- counts %>% distinct(condition, bin, condition_key) %>% 
  left_join(coords, by = 'bin') %>%
  mutate(sample_id = condition, 
         grid_key = condition_key) %>%
  arrange(grid_key)

#============================================================
stopifnot('binned counts does not match with metadata' = all(rownames(prop) == bin_meta$grid_key))
#============================================================
# Verify means across all bins:
tapply(prop[, 'cMESHyp'], bin_meta$condition, mean)
tapply(prop[, 'Vasc'],    bin_meta$condition, mean)
tapply(prop[, 'cOXPHOS'], bin_meta$condition, mean)
#============================================================
prop_centered <- scale(prop, center = TRUE, scale = TRUE)
#============================================================
par(mfrow =c(2,2))
boxplot(prop_centered[, 'cOXPHOS']~ bin_meta$condition, main = 'cOxphos')
boxplot(prop_centered[, 'cMESHyp']~ bin_meta$condition, main = 'cMesHyp')
boxplot(prop_centered[, 'Vasc']~ bin_meta$condition, main = 'Vasc')
boxplot(prop_centered[, 'cMESInf']~ bin_meta$condition, main = 'cMesInf')
#============================================================
# filter at the bin level
# Require a bin to have sufficient spots in BOTH conditions (or relax threshold for the minor condition)
# valid_bins <- counts %>%
#   distinct(condition_key, total_cells) %>%
#   tidyr::separate(condition_key, into = c("condition", "bin"), sep = "___") %>%
#   tidyr::pivot_wider(names_from = condition, values_from = total_cells, values_fill = 0) %>%
#   # Example: >= 80 total spots, and at least 20 spots in each condition so proportions are well-defined
#   filter((GBM + IDHmutant) >= 80, GBM >= 20, IDHmutant >= 20) %>%
#   pull(bin)
# 
# bin_meta <- bin_meta %>% filter(bin %in% valid_bins)
# bin_meta$bin <- droplevels(bin_meta$bin)
(N <- length(levels(bin_meta$bin)))
#============================================================
# # Variance-stabilizing transformation for multinomial proportions, centered per state
# prop <- prop[as.character(bin_meta$condition_key), ]
# prop_vst <- asin(sqrt(prop))
# prop_centered <- scale(prop_vst, center = TRUE, scale = FALSE)



#============================================================
all(rownames(prop_centered) == bin_meta$condition_key)
#========================================================================
## build the full unobserved regular spatial grid (n_slides x N)
#========================================================================
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
image(projection_matrix)
dim(projection_matrix)
#=========================================
## prepare the coregionalization matrix
#=========================================
W <- model.matrix(~0 + levels(bin_meta$condition))
B <- tcrossprod(W) 
image(B)

noise_diag <- c(1, 2)
## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + bin_meta$sample_id)
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)
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
data_env$donor_identity <- cov_donor
data_env$noise_diag <- noise_diag # mapping from individual samples to condition
data_env$KERNEL = 'RBF'
data_env$FIT_DONOR <- FALSE
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
data_env$RETURN_FITS <- TRUE
data_env$kernel_params <- NULL # if my parameters are already fit, I don't want to fit likelihood again

## save variables to memory
write_rds(data_env, file = 'output/glioma_condition_kmeans_data_env.rds')
#============================================================
## source implementation modules
#============================================================
res.ARGO <- fit_ARGO(data_env)
res.ARGO$params
write_rds(res.ARGO, file = 'output/WT_vs_mutant_kmeans_results.rds')
#============================================================
# data_env <- read_rds('output/glioma_condition_data_env.rds')
# list2env(as.list(data_env), envir = globalenv())
# res.ARGO <- read_rds('output/WT_vs_mutant_condition_results.rds')
#=================================================
## visualize results
#=================================================
funGP_fits <- res.ARGO$resObjects
funGP_res <- funGP_fits %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation)) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 
funGP_res
#=================================================
## where the significant deviation lies
#=================================================
ggplot(funGP_res %>% arrange(desc(stats)), 
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
    aes(label = gene),
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


v_thresh <- 1
h_thresh <- -log(0.05)
funGP_res %>% mutate(status = case_when(
  neglogp > h_thresh & deviation > v_thresh ~ 'Upreg',
  neglogp > h_thresh & deviation < -v_thresh ~ 'Downreg',
  TRUE ~ 'Not Significant'
)) %>% arrange(status == 'Not Significant') %>%
  ggplot(aes(x = deviation, y = neglogp)) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status), alpha = 0.35, size = 2, shape = 20) +
  scale_color_manual(
    values = c(
      "Upreg" = "#e41a1c",   # Professional Red
      "Downreg" = "#377eb8",  # Professional Blue
      "Not Significant" = "gray60"     # Muted Gray
    )
  ) + 
  geom_text_repel(
    data = . %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 10))},
    aes(label = gene),
    size = 5, 
    family = 'Arial',
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray", 
    max.overlaps = 30
  ) +
  # 2. Horizontal significance threshold line: -log10(0.05) ~ 1.301
  geom_hline(yintercept = h_thresh, color = "red", linetype = "dashed", linewidth = .75) +
  # 3. Vertical effect size threshold lines
  geom_vline(xintercept = c(-v_thresh, v_thresh), color = "red", linetype = "dashed", size = .75) +
  # 4. Axis Labels
  labs(
    x = "max deviation (LFC)",
    y = expression(-log~"analytic")
  ) +
  # xlim(c(-3.5, 3.5)) + 
  # 5. Clean, professional theme minimalization
  theme_bw(base_family = 'Arial') +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size=14)
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1))
  )
