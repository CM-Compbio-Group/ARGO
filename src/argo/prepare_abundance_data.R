#' prepare cell type abundance data for modeling
source('src/argo/00-load-packages.R')
CHECK_SPARSITY <- FALSE

ad <- anndata::read_h5ad('~/Desktop/integrated_anndata.h5ad')

## only keep the columns of interest
obs <-  ad$obs %>% rownames_to_column('cell_id') %>% 
  as_tibble() %>% select(cell_id,  donor_id, clust_annot,
                         condition = age, sample_id = batch,
                         pseudotime = binned_isodepth) %>% 
  mutate(condition = factor(condition, levels = c('4wk', '24wk', '90wk')),
         sample_id = factor(sample_id),
         donor_id = factor(donor_id, levels = gtools::mixedsort(unique(donor_id))))

#============================================================
## find what resolution to fit data at by looking at gene expression sparsity levels
#============================================================

find_sparsity <- function(N, obs){
  abundance <- obs %>% mutate(bin = ntile(pseudotime, N)) %>%
    mutate(grid_key = paste(sample_id, bin, sep = '___')) %>%
    count(sample_id, bin, grid_key, clust_annot, name = 'n_cells') %>%
    group_by(sample_id, bin) %>% 
    mutate(total_cells = sum(n_cells),
           proportion = n_cells / total_cells) %>% 
    ungroup()
  
  prop <- abundance %>% select(grid_key, clust_annot, proportion) %>%
    pivot_wider(
      names_from = clust_annot,           # Cell types → column names
      values_from = proportion,        # Proportions → cell values
      values_fill = 0,                 # NA → 0 (missing cell types)
      names_sort = TRUE                # Alphabetize cell type columns
    ) %>% column_to_rownames('grid_key') %>% as.matrix()
  
  return(rowSums(prop == 0))
}

if(CHECK_SPARSITY){
  n_range <- seq(10,90,10)
  sparsity <- lapply(n_range, find_sparsity, obs)
  sapply(sparsity, function(x) sum(x > 0))
  
  par(mfrow = c(3,3))
  for(i in seq_along(n_range)){
    hist(sparsity[[i]],  xlab = 'number of sparse cell-types per bin', 
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

obs <- obs %>% 
  mutate(bin = factor(ntile(pseudotime, n = N))) %>% 
  group_by(bin) %>% mutate(coord_x = mean(pseudotime)) %>% 
  ungroup() %>%
  mutate(grid_key = paste(sample_id, bin, sep = '___'), 
         condition_key = paste(condition, bin, sep = '___')) %>%
  mutate(grid_key = factor(grid_key, 
                           levels = gtools::mixedsort(unique(grid_key))),
         condition_key = factor(condition_key, 
                                levels = gtools::mixedsort(unique(condition_key))))


# cell type abundance data
abundance <- obs %>% 
  dplyr::count(grid_key, clust_annot, name = 'n_cells') %>%
  group_by(grid_key) %>%
  mutate(total_cells = sum(n_cells),
         proportion = n_cells / total_cells) %>%
  ungroup()

prop <- abundance %>% select(grid_key, clust_annot, proportion) %>%
  pivot_wider(
    names_from = clust_annot,           # Cell types → column names
    values_from = proportion,        # Proportions → cell values
    values_fill = 0,                 # NA → 0 (missing cell types)
    names_sort = TRUE                # Alphabetize cell type columns
  ) %>% column_to_rownames('grid_key') %>% as.matrix()

# pheatmap(t(prop), cluster_cols = T)
#=============================
# prepare bin level metadata
#=============================
bin_meta <- obs %>% select(-c(cell_id, pseudotime, clust_annot)) %>% 
  distinct() %>% arrange(grid_key)

stopifnot('binned counts does not match with metadata' = all(rownames(prop) == bin_meta$grid_key))
#===============================================
# trends for olig state proportions across bins
#===============================================
## convert to mean scaled counts in each condition x bin 
n_counts <- bin_meta %>% dplyr::count(condition_key)
mean_abx <- rowsum(prop, bin_meta$condition_key) 
## divide by frequency of condition key in bin_meta
all(rownames(mean_abx) == n_counts$condition_key)
mean_abx <- mean_abx / n_counts$n 
mean_abx <- mean_abx %>% 
  as.data.frame() %>% 
  rownames_to_column('condition_key') 

mean_abx <- pivot_longer(mean_abx, cols = -condition_key, 
                         names_to = 'cell_type', 
                         values_to = 'counts') %>% 
  left_join(bin_meta %>% select(condition, coord_x, condition_key),
            by = 'condition_key', multiple = 'first') %>% 
  filter( cell_type %in% c("Olig-1", "Olig-2", "Olig-3")) %>% 
  group_by(condition_key) %>%
  mutate(total_frac = sum(counts)) %>%
  ungroup()

# pal43 <- colorRampPalette(RColorBrewer::brewer.pal(6, "Set2"))(3)
# pal43 <- c('gray85', 'slategray', 'gray25')
pal43 <- c('lightpink', 'palevioletred', 'violetred','darkred')

p <- ggplot(mean_abx, aes(x = coord_x, 
                          y = counts,  # Or n_cells
                          color = cell_type, 
                          group = cell_type)) +
  geom_line(linewidth = .5, alpha = 0.75, linetype=2) +
  geom_line(aes(x = coord_x, y = total_frac, 
                color = "Total Fraction", group = 'condition_key'
                ), 
            linewidth = .5, alpha = 0.95, linetype = 2) +
  facet_wrap(~ interaction(condition), 
             scales = "free_y", ncol=1) +
  labs(
    # title = "Mean Cell Type countsortions vs. Isodepth Bins by Age",
    x = "Isodepth", y = "Scaled Abundance",
    color = "Cell Type") +
  guides(color = guide_legend(nrow=2,
    override.aes = list(
      linewidth = rep(.75, 4), # Thicker for Total Fraction (last item)
      linetype  = c(2, 2, 2, 2)             # Solid (1) for Total Fraction, Dashed (2) for others
    )
  )) +
  scale_color_manual(values = pal43) +
  theme_classic(base_family = 'Arial', 
                base_size = 8) +
  # guides(color = guide_legend(ncol = 3)) +  # Set number of columns
  theme(
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom", 
        plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        
        strip.text = element_text(size = 8, face = "plain"), 
        # strip.background = element_blank(),           # Optional: Cleans up facet header boxes
        
        # Legend text and title formatting
        legend.title = element_text(size = 8, face = "plain"),
        legend.text = element_text(size = 8, face = "plain"), 
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0), # Pulls legend closer to plot
        legend.key.size = unit(3, "mm")                 # Shrinks the color boxes
  ) + scale_y_continuous(
    expand = expansion(mult = c(0.025, 0.15))) +
  coord_cartesian(clip = "off") 

print(p)

ggsave(plot = p, filename = glue('figs/Olig_fractions.pdf'), 
       width = 3.5, height =4, units = 'in')

# write_rds(mean_abx, file = 'output/mean_olig_abundance.rds') # save to make the abundance polots 
#====================================
# filter out bins with high sparsity
#====================================
bin_meta <- bin_meta %>% 
  mutate(sparsity = rowSums(prop == 0)) 

## see if sparsity is related to bin identity
filter(bin_meta, sparsity > 25) %>% with(table(bin))
# data is not missing at random; old mice have fewer cell type diversity
# this could just be true biological signal too. 

# i'll remove outlier bins with over 25 missing cell types
sum(bin_meta$sparsity > 25) / nrow(bin_meta) # about 4% of bins
bin_meta <- filter(bin_meta, sparsity < 25)
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
stopifnot('counts matrix is not aligned with metadata' = all(rownames(prop) == bin_meta$grid_key))
#============================================================
# center the abundance ata
#============================================================
prop <- scale(prop, center = T, scale = T)
# pheatmap(prop)

plot(prop[, 'Olig-1']~ bin_meta$condition_key)
plot(prop[, 'Olig-2']~ bin_meta$condition_key)
plot(prop[, 'Olig-3']~ bin_meta$condition_key)
plot(prop[, 'Endo-1']~ bin_meta$condition_key)
plot(prop[, 'Endo-3']~ bin_meta$condition_key)
plot(prop[, 'Micro-3']~ bin_meta$condition_key)
#============================================================

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
data_env$cov_donor <- cov_donor
data_env$noise_diag <- noise_diag
data_env$KERNEL = 'RBF'
data_env$FIT_DONOR <- FALSE
data_env$FIT_CONDITIONAL_NOISE <- TRUE
data_env$LAMBDA <- 1 # regularization for noise terms in GP fits

data_env$CONTRAST_VECTOR <- c(-1, 1, 0)
data_env$CONTRAST <- 'adult_vs_juvenile'
data_env$RETURN_CONF_BANDS = TRUE # compute conf bands only if I am interested in recovering intervals
data_env$conf_level_alpha = 0.95
data_env$n_monte_carlo_simulations = 1000 # no.of monte carlo samples
data_env$features_to_fit <- colnames(prop)
data_env$test_coords <- NULL
data_env$PLOT_POSTERIOR <- FALSE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$RETURN_FITS <- TRUE
# data_env$kernel_params <- NULL
data_env$kernel_params <- res.ARGO$params # if my parameters are already fit, I don't want to fit likelihood again

## save variables to memory
write_rds(data_env, file = 'output/abundance_data_env.rds')
#============================================================
## source implementation modules
#============================================================
res.ARGO <- fit_ARGO(data_env)
write_rds(res.ARGO, file = glue('output/{data_env$CONTRAST}_diff_abundance_results.rds'))
#============================================================
## visualize results
#============================================================
data_env$CONTRAST <- 'old_vs_adult'
res.ARGO <- read_rds(glue('output/old_vs_adult_diff_abundance_results.rds'))

## compute exceedence metrics
exc.metrics <- do.call('rbind', lapply(res.ARGO$resObjects, compute_exceedence_metrics, isodepth_))
## aggregate results from functional testing
funGP_res <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  bind_cols(exc.metrics) %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         net_auec = auec_up - auec_down,
         net_poe = poe_up - poe_down) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(stats)) %>% 
  mutate(Rank = row_number()) 

#============================================================
# Trumpet plot
pdf(glue('figs/{data_env$CONTRAST}_abundance_ranked_plot.pdf'), width = 6, height = 6)
plot(ggplot(funGP_res %>% arrange(desc(stats)), aes(x = Rank, y = stats)) +
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
    title = glue("{data_env$CONTRAST}")
  ) +
  theme_classic(base_family = 'Arial') +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size=14)
  ))
dev.off()

#============================================================
# volcano plot
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
    data = . %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))},
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
  geom_hline(yintercept = h_thresh, color = "red", linetype = "dashed", size = .75) +
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
