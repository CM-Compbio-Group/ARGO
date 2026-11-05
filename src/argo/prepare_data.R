#' prepare the brain aging data for fitting with ARGO
source('src/argo/00-load-packages.R')

CHECK_SPARSITY <- FALSE

ad <- anndata::read_h5ad('~/Desktop/integrated_anndata.h5ad')
counts <- as.matrix(ad$layers['counts'])

## only keep the columns of interest
metadata <-  ad$obs %>% rownames_to_column('cell_id') %>% 
                      as_tibble() %>% select(cell_id,  donor_id, 
                      condition = age, sample_id = batch,
                      pseudotime = binned_isodepth) %>% 
  mutate(condition = factor(condition, levels = c('4wk', '24wk', '90wk')),
         donor_id = factor(donor_id, levels = gtools::mixedsort(unique(donor_id))))


stopifnot('counts ordering does not match with metadata' = all(rownames(counts) == metadata$cell_id))

#======================================
if(CHECK_SPARSITY){
  source('src/R/check_sparsity.R')
}

#============================================================
## fitting data with 30 bins, I will marginalize over bins that have high sparsity
#============================================================
N <- 30 # number of bins to fit data

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
dim(projection_matrix)
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
data_env$cov_donor <- cov_donor
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
data_env$features_to_fit <- colnames(bin_cts)
data_env$test_coords <- NULL
data_env$kernel_params <- NULL
data_env$GET_POSTERIOR <- FALSE
data_env$PLOT_POSTERIOR <- FALSE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$RETURN_FITS <- TRUE

write_rds(data_env, file = 'output/expression_data_env.rds')
# res.ARGO <- fit_ARGO(data_env)


#============================================================
## timed run
#============================================================

data_env <- read_rds('output/expression_data_env.rds')
data_env$GET_POSTERIOR <- FALSE
data_env$RETURN_FITS <- TRUE
data_env$FIT_DONOR<-TRUE
data_env$features_to_fit <- colnames(data_env$counts)

start_time <- Sys.time()
res.ARGO <- fit_ARGO(data_env)
end_time <- Sys.time()
difftime(end_time, start_time)

print(glue('Time taken to run with N={data_env$N} is:{difftime(end_time, start_time, units = "secs")}'))
write_rds(res.ARGO, file = glue('output/{data_env$CONTRAST}_N30_w_donor_diff_testing.rds'))

#=====================================================
funGP.df <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())

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

