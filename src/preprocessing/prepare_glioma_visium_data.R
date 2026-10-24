#' AddModuleScore for ivy niche enrichment scores
#' add a soft-plus transform, rather than clipping scores at zero, to retain 
#' identifiability across bins
#' scale the inputs before generating bins
#' get pooled abundance data by condition 
#' filter out bins with very few spots -- thresholding at 100 spots here.
#' define the projection and donor_covariance matrices.
#' scale abundance data before fitting. 

source('src/argo/00-load-packages.R')
library(Seurat)
library(RColorBrewer)


brain <- read_rds('~/Desktop/brain-seurat-objects/visium.rds')
source('src/R/bin_kdtree_adaptive.R')

ivy <- levels(brain$ivy_niche)

# xx <- brain@meta.data %>% dplyr::select(all_of(ivy))
## clamp the enrichment scores at zero
# xx <- pmax(as.matrix(brain@meta.data %>% dplyr::select(all_of(ivy))), 0)

## convert to relative scores for spots
raw_scores <- as.matrix(brain@meta.data[, ivy])
# Temperature tau = 0.15 controls sharpness toward the 4 vertices
beta <- 0.05
scores_softplus <- scale(beta * log1p(exp(raw_scores/beta)))

par(mfrow = c(2, 2))
for (col in colnames(scores_softplus)) {
  boxplot(scores_softplus[, col] ~ brain$ivy_niche, 
          main = paste(col, "(Softplus)"),
          xlab = "Assigned Niche",
          ylab = col,
          outline = FALSE)
}

par(mfrow = c(1,1))

## bin the data 
brain$bins <- bin_kdtree_adaptive(scores_softplus, target_spot_count = 1000)
(N <- length(levels(brain$bins)))

df <- scores_softplus %>% 
  bind_cols(brain@meta.data %>% dplyr::select(
                                        ivy_niche, 
                                        MP = annotation1,
                                        condition = subtype, 
                                        sample_id = orig.ident,
                                        bin = bins))  %>% 
                         mutate(grid_key = paste(sample_id, bin, sep = '___'), 
                                condition_key = paste(condition, bin, sep = '___'),
                                condition = factor(condition)) %>%
                         mutate(sample_id = factor(sample_id, 
                                            levels = gtools::mixedsort(unique(sample_id))),
                                bin = factor(bin, levels = sort(unique(bin))), 
                                grid_key = factor(grid_key, 
                                            levels = gtools::mixedsort(unique(grid_key))),
                                condition_key = factor(condition_key, 
                                                  levels = gtools::mixedsort(unique(condition_key)))) %>%
  arrange(grid_key)  %>% 
  rownames_to_column('spot_id') %>% 
  as_tibble()
                        

bin_cols <- colorRampPalette(brewer.pal(9, "Set3"))(N)
plotly::plot_ly(data = df, 
                x = ~CT, 
                y = ~CTpan, 
                z = ~LE, 
                color = ~condition_key, 
                colors = bin_cols,
                type = "scatter3d", 
                marker = list(size = 2.5, alpha=.5),
                mode = "markers")


## filter out bins that are extremely skewed b/w conditions
bin_tab <- table(df$bin, df$condition)
prop_bins <- bin_tab[, "GBM"] / rowSums(bin_tab)

# filter out bins with extreme skew b/w conditions
# keep_bins <- which((prop_bins > .05) & (prop_bins < .95))
# df <- filter(df, bin %in% keep_bins)
# df$bin <- droplevels(df$bin)
# (N <- length(levels(df$bin)))
#============================================================
## prepare environment for running argo
#============================================================
# I am running diff tests b/w conditions; pooling across all slides per condition
# have no use for grid_key here -- all information is summarized at level of condition_key
# coords per bin per condition; counts per bin per condition
counts <- df %>%
  dplyr::count(condition_key, 
              MP, name = 'n_cells') %>% 
  group_by(condition_key) %>%
  mutate(total_cells = sum(n_cells),
         proportion = n_cells / total_cells) %>%
  ungroup()

prop <- counts %>% select(condition_key, MP, proportion) %>%
  arrange(condition_key) %>% 
  pivot_wider(
    names_from = MP,           # Cell types → column names
    values_from = proportion,        # Proportions → cell values
    values_fill = 0,                 # NA → 0 (missing cell types)
    names_sort = TRUE                # Alphabetize cell type columns
  ) %>% 
  column_to_rownames('condition_key') %>%
  as.matrix()

# pheatmap(prop, cluster_rows = T)

# each bin has a coordinate in space
coords <- df %>% 
  group_by(bin) %>% 
  summarise(across(all_of(ivy), \(x) mean(x, na.rm = T)))

coords

bin_meta <- left_join(left_join(coords, df %>% 
  distinct(condition, bin, condition_key), 
  by = 'bin'), counts %>% distinct(condition_key, total_cells), 
  by = 'condition_key') %>% arrange(condition_key)

stopifnot('binned counts does not match with metadata' = all(rownames(prop) == bin_meta$condition_key))

## now filter out bins with very few cells
# hist(bin_meta$total_cells)
# quantile(bin_meta$total_cells, c(.1,.2, .25,.5))
# sum(bin_meta$total_cells < 5)

# remove binxslide combinations with fewer than 10 cells
sum(bin_meta$total_cells < 100)

# drop bins with fewer than 100 spots
bin_meta <- bin_meta %>% filter(total_cells > 100)
prop <- prop[bin_meta$condition_key, ]
#========================================================================
## build the full unobserved regular spatial grid (n_slides x N)
#========================================================================
# sample-1 is condition-1; sample-2 is condition-2
bin_meta <- bin_meta %>% mutate(sample_id = fct_recode(condition, 
                                                       'sample_1' = 'GBM',
                                                       'sample_2' = 'IDHmutant'), 
                                grid_key =  condition_key)

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
image(projection_matrix)
#=========================================
## prepare the coregionalization matrix
#=========================================
W <- model.matrix(~0 + levels(bin_meta$sample_id))
B <- tcrossprod(W) 
# image(B)

noise_diag <- c(1, 2)
## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + bin_meta$sample_id)
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)

stopifnot('counts matrix is not aligned with metadata' = all(rownames(prop) == bin_meta$condition_key))
#============================================================
# center the abundance ata
#============================================================
prop <- scale(prop, center = T, scale = T)


# as_tibble(prop) %>% 
#   mutate(condition = bin_meta$condition) %>% 
#   group_by(condition) %>% 
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

par(mfrow = c(1,2))
plot(prop[, 'Vasc']~ bin_meta$condition)
plot(prop[, 'Neuron']~ bin_meta$condition)

plot(prop[, 'Vasc']~ bin_meta$total_cells, pch = 20, col = factor(bin_meta$condition))
plot(prop[, 'Neuron']~ bin_meta$total_cells, pch = 20, col = factor(bin_meta$condition))

plot(prop[, 'Neuron']~ bin_meta$LE, pch = 20, col = factor(bin_meta$condition))
plot(prop[, 'Vasc']~ bin_meta$CTmvp, pch = 20, col = factor(bin_meta$condition))


plot(prop[, 'cMESHyp']~ bin_meta$condition)
plot(prop[, 'cMESInf']~ bin_meta$condition)
plot(prop[, 'cMESVasc']~ bin_meta$condition)
plot(prop[, 'cProlifMetab']~ bin_meta$condition)
plot(prop[, 'cOXPHOS']~ bin_meta$condition)
plot(prop[, 'Oligo']~ bin_meta$condition)
plot(prop[, 'cNPC']~ bin_meta$condition)
plot(prop[, 'cOPC1']~ bin_meta$condition)
par(mfrow=c(1,1))

#============================================================
#the coordinate axis for input
isodepth_ <- bin_meta %>% select(all_of(ivy)) %>% distinct() %>% as.matrix()
# bin_meta$grid_key <- bin_meta$condition_key # need grid key for core code
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
data_env$noise_diag <- noise_diag # mapping from individual samples to condition
data_env$KERNEL = 'RBF'
data_env$FIT_DONOR <- FALSE
data_env$FIT_CONDITIONAL_NOISE <- TRUE
data_env$LAMBDA <- 1 # regularization for noise terms in GP fits

data_env$CONTRAST_VECTOR <- c(1, -1)
data_env$CONTRAST <- 'WT_vs_Mutant'
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
write_rds(data_env, file = 'output/glioma_condition_data_env.rds')
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
source('src/R/compute_exceedence_metrics.R')
source('src/R/plot_conf_bands.R')

res.ARGO <- fit_ARGO(data_env)
write_rds(res.ARGO, file = 'output/WT_vs_mutant_condition_results.rds')

# I don't need to fit the model within each slide separately.
# I can just jointly model across all slides and split by condition -- not modeling donor effects

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


v_thresh <- .01
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
    data = . %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))},
    aes(label = gene),
    size = 5, 
    family = 'Arial',
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray", 
    max.overlaps = 20
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
