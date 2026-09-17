#' **0B-load-data-for-modeling** : load the bin-level counts data-metadata object and prepares all variables necessary for ICM modeling with funGP
#' organizes the data in proper form, and derives the supplemental variables like co-regionalization matrix and noise matrices necessary for fitting GP models
#' input : processed bin-level counts data -- counts are centered across all observations, *not scaled*
#'         data optionally has reflective padding, if desired


# read processed bin-level counts data
bin_meta <- read_rds('data/processed/bin_meta.rds')
#--------------------------------
#--------------------------------
# extract variable names
genes <- colnames(bin_meta)[9:382]
## extract slides in strict condition--donor order
slide_df <- bin_meta %>% 
  distinct(age, donor_id, batch) %>%
  arrange(age, donor_id, batch)

ordered_slides <- slide_df$batch
slides_per_age <- slide_df %>% dplyr::count(age) %>% pull(n)
n_slides <- nrow(slide_df)

## load the coordinate points (test data points)
## observed data must include the padding
x_observed <- sort(unique(bin_meta$coord_x)); N <- length(x_observed)
load('data/temp/test_domain.rdata')
#--------------------------------------------------------------------
## build the full unobserved regular spatial grid (n_slides x N)
#--------------------------------------------------------------------
ordered_bins <- levels(bin_meta$bin)

full_grid <- expand_grid(
  batch = ordered_slides, ## first arguement varies slowest
  coord_x = ordered_bins,
  ) %>% left_join(slide_df, by = 'batch') %>%
  mutate(grid_key = paste(batch, coord_x, sep = '___'))

## make the same variable with bin_meta
bin_meta <- bin_meta %>% mutate(grid_key = paste(batch, bin, sep = '___')) %>%
  arrange(match(grid_key, full_grid$grid_key))
  
## construct the permutation matrix P_obs
matched_indices <- match(bin_meta$grid_key, full_grid$grid_key)
n_obs <- nrow(bin_meta); n_full <- nrow(full_grid)
P_obs <- matrix(0, nrow = n_obs, ncol = n_full)
for(i in seq_len(n_obs)){
  P_obs[i, matched_indices[i]] <- 1
}
#--------------------------------------------------------------------
## construct the co-regionalization and donor-covariance matrix
#--------------------------------------------------------------------
W <- model.matrix(~0 + factor(slide_df$age))
rho <- 1
B <- rho * tcrossprod(W) + (1-rho)*diag(n_slides)
image(B)

noise_diag <- as.numeric(factor(slide_df$age))

## donor-intercept structure in observed space
Z_donor <- model.matrix(~0 + factor(bin_meta$donor_id))
cov_donor <- tcrossprod(Z_donor)
image(cov_donor)
#-------------------------------------------------------------
## visualize expression of gene from aligned data
#-------------------------------------------------------------
# plot(bin_meta$C4b~ bin_meta$age)
#------------------------------------------------
## set parameters and values for downstream plotting
#------------------------------------------------
## define cov kernel exponents
## pairwise distances and dist-squared
# r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
# r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
# r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2
#------------------------------------------------
## pick genes to test fits on
#------------------------------------------------
# sig_genes <- vroom('data/processed/sig_de_genes.csv')
# genes_to_plot <- filter(sig_genes, shape_cluster == 4) %>% pull(gene)
