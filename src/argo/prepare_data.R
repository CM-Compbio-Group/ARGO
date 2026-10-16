#' prepare the brain aging data for fitting with ARGO

library(tidyverse)
library(anndata)
library(glue)
library(patchwork)
library(ComplexHeatmap)
library(ggrepel)

CHECK_SPARSITY <- FALSE

ad <- read_h5ad('~/Desktop/integrated_anndata.h5ad')
counts <- as.matrix(ad$layers['counts'])

## only keep the columns of interest
metadata <-  ad$obs %>% rownames_to_column('cell_id') %>% 
                      as_tibble() %>% select(cell_id,  donor_id, 
                      condition = age, sample_id = batch,
                      pseudotime = binned_isodepth) %>% 
  mutate(condition = factor(condition, levels = c('4wk', '24wk', '90wk')),
         donor_id = factor(donor_id, levels = gtools::mixedsort(unique(donor_id))))


stopifnot('counts ordering does not match with metadata' = all(rownames(counts) == metadata$cell_id))
#============================================================
## find what resolution to fit data at by looking at gene expression sparsity levels
#============================================================

find_sparsity <- function(N, counts, metadata){
  bins <- ntile(metadata$pseudotime, N)
  slide_bin <- factor(paste(metadata$sample_id, bins, sep = '_'))
  bin_cts <- rowsum(counts, group = slide_bin)
  return(rowSums(bin_cts == 0))
}

if(CHECK_SPARSITY){
  n_range <- seq(10,90,10)
  sparsity <- lapply(n_range, find_sparsity, counts, metadata)
  sapply(sparsity, function(x) sum(x > 0))
  
  par(mfrow = c(3,3))
  for(i in seq_along(n_range)){
    hist(sparsity[[i]],  xlab = 'number of sparse genes per bin', 
         main = glue('Sparsity with {n_range[i]} isodepth bins'))
  }
  
  sparsity_quantiles <- lapply(sparsity, function(x) quantile(x, c(.5, .75, .9, .99)))
  
  t <- as.matrix(do.call('rbind', sparsity_quantiles))
  rownames(t) <- n_range
  print(t)
}




#============================================================
## fitting data with 30 bins, I will marginalize over bins that have high sparsity
#============================================================
N <- 30 # number of bins to fit data
M <- N # testing at same resolution

metadata <- metadata %>% 
  mutate(bin = factor(ntile(pseudotime, n = N))) %>% 
  group_by(bin) %>% mutate(coord_x = mean(pseudotime)) %>% 
  ungroup() %>%
  mutate(grid_key = paste(sample_id, bin, sep = '___'), 
         condition_key = paste(condition, bin, sep = '___')) %>%
  mutate(grid_key = factor(grid_key, 
                             levels = gtools::mixedsort(unique(grid_key))),
         condition_key = factor(condition_key, 
                             levels = gtools::mixedsort(unique(condition_key))))


# prepare the pseudobulk counts data
bin_cts <- rowsum(as.matrix(counts), group = metadata$grid_key)

#=============================
# prepare bin level metadata
#=============================
bin_meta <- metadata %>% select(-c(cell_id, pseudotime)) %>% 
  distinct() %>% arrange(grid_key)

stopifnot('binned counts does not match with metadata' = all(rownames(bin_cts) == bin_meta$grid_key))
#====================================
# filter out bins with high sparsity
#====================================
bin_meta <- bin_meta %>% 
            mutate(sparsity = rowSums(bin_cts == 0)) %>%
            dplyr::filter(sparsity < 3)

bin_cts <- bin_cts[bin_meta$grid_key, ]

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

#=========================================
## prepare the coregionalization matrix
#=========================================
slide_df <- bin_meta %>% 
            distinct(sample_id, donor_id, condition)

W <- model.matrix(~0 + slide_df$condition)
B <- tcrossprod(W) 
image(B)

noise_diag <- as.numeric(slide_df %>% 
                           pull(condition))

## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + bin_meta$donor_id)
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)


isodepth_ <- unique(bin_meta$coord_x)
stopifnot('counts matrix is not aligned with metadata' = all(rownames(bin_cts) == bin_meta$grid_key))

#============================================================
# normalize and center the counts data
#============================================================
bin_cts <- scale(t(DESeq2::varianceStabilizingTransformation(t(bin_cts))), 
                 center = T, scale = F)

#============================================================
# create environment to pass these variables for fitting
#============================================================
data_env <- new.env()
data_env$counts <- bin_cts
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
data_env$CONTRAST_VECTOR <- c(0, -1, 1)
data_env$CONTRAST <- 'old_vs_adult'
data_env$RETURN_CONF_BANDS = TRUE # compute conf bands only if I am interested in recovering intervals
data_env$conf_level_alpha = 0.95
data_env$n_monte_carlo_simulations = 1000 # no.of monte carlo samples
# data_env$features_to_fit <- colnames(bin_cts)
data_env$features_to_fit <- top_diff
data_env$test_coords <- NULL
data_env$kernel_params <- res.ARGO$params[top_diff, ]


#============================================================
## read downstream functions
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

funGP.df <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())

write_rds(res.ARGO, file = 'output/old_vs_adult_N30_diff_testing_results.rds')

#============================================================
## make dotplots for the top differential genes
#============================================================
sce <- SingleCellExperiment(assays = list(counts = t(counts)), 
                            colData = metadata)
sce <- scater::logNormCounts(sce)

dot_df <- as.data.frame(colData(sce)[, c("condition", "donor_id", "sample_id")]) %>%
  mutate(sample_id = paste('slide', sample_id, sep = '_')) %>%
  bind_cols(as.data.frame(t(as.matrix(assay(sce, "logcounts")[genes_to_plot, ]) ))) %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  group_by(gene, condition, donor_id, sample_id) %>%
  summarize(
    pct_exp = mean(expr > 0) * 100,
    mean_exp = mean(expr),
    .groups = "drop"
  ) %>%
  # Z-score mean expression per gene so all genes share a comparable color bar
  group_by(gene) %>%
  mutate(scaled_mean = as.numeric(scale(mean_exp))) %>%
  ungroup()

p <- ggplot(dot_df, aes(x = sample_id, y = factor(gene, levels = rev(genes_to_plot)))) +
  geom_point(aes(size = pct_exp, color = scaled_mean)) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(.25, 5), name = "% Expressing") +
  scale_color_gradient2(low = 'blue3', high='red3', mid = 'lightgrey', 
                        name = "Z-Scored\nMean Expr")+
  labs(x = "", y = "Gene") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

ggsave(plot=p, filename = 'figs/top_diff_dotplot.pdf', width = 10, height = 6)

