# a function for computing Jaccard index
union_length <- function(l, r) {
  if (length(l) == 0) return(0)
  
  valid <- l < r
  l <- l[valid]
  r <- r[valid]
  
  if (length(l) == 0) return(0)
  
  o <- order(l)
  l <- l[o]
  r <- r[o]
  
  sum(pmax(0, r - pmax(l, c(-Inf, head(cummax(r), -1)))))
}

compute_jaccard_sim <- function(fdr, cell_type_range,
                                true_l, true_r, fdr_threshold) {
  sig <- !is.na(fdr) & fdr < fdr_threshold
  
  l <- cell_type_range$min_pseudotime[sig]
  r <- cell_type_range$max_pseudotime[sig]
  
  # Total length of predicted intervals, excluding gaps
  pred_len <- union_length(l, r)
  
  # Total overlap with true interval
  intersect_len <- union_length(
    pmax(l, true_l),
    pmin(r, true_r)
  )
  
  # Inclusion-exclusion
  union_len <- pred_len + (true_r - true_l) - intersect_len
  
  if (union_len == 0) return(NA_real_)
  intersect_len / union_len
}

#======================================================
#======================================================
# read simulated data for this experiment iteration and unload the list into variables
zeallot::`%<-%`(c(sim_counts, sim_meta, genes, effect_info), 
                read_rds(glue('testing/data/temp/{dataset}_w_replicates_interval_perturbation.rds')))

## update effect_info to track 20% FC as the signal thresholds
# effect_info <- effect_info %>% 
#   dplyr::rename(interval_left = interval_left_01,
#                 interval_right = interval_right_01)
# 
# effect_info <- left_join(effect_info, rel_length, by = 'bandwidth') %>% 
#   mutate(across(rel_length, ~ round(., 3))) %>% 
#   dplyr::rename(`relative length` = rel_length,
#                 `signal strength` = amplitude)
# 
# # get metadata for different scales of data representation
zeallot::`%<-%`(c(cell_level_meta, bin_level_meta,
                  sample_level_meta, cell_type_meta), 
                prepare_metadata(sim_meta, N))

cell_type_range <- cell_level_meta %>% group_by(cell_type) %>%
  summarise(
    min_pseudotime = quantile(pseudotime, 0.02, na.rm = TRUE),
    max_pseudotime = quantile(pseudotime, 0.98, na.rm = TRUE),
    .groups = "drop"
  )

## map cell types to pseudotime intervals
cluster_intervals <- as_tibble(sce@colData) %>%  
  mutate(cell_type = str_replace_all(cell_type, 
                     c(' ' = '_', '-' = '_'))) %>% 
  group_by(cell_type) %>%
  summarize(left_ = quantile(pseudotime, .02),
            right_ = quantile(pseudotime, .98))

## convert to pseudobulk counts
celltype_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_cluster)
## assert that sample IDs match
stopifnot(all(rownames(celltype_counts) == cell_type_meta$sample_cluster))
clusters <- unique(cell_type_meta$cell_type)
#======================================================
## perform differential testing
#======================================================
diff.testing.methods <- list(
  DESeq2 = run_deseq, 
  edgeR = run_edger, 
  Limma = run_limma
)
## capture the current environment to pass on to differential testing functions
env <- environment() 

## call differential testing passing the local environment
p.vals <- map(diff.testing.methods, ~ .x(env))
fdr.list <- lapply(p.vals, function(res) apply(res, 2, p.adjust, method = 'BH'))

write_rds(fdr.list, 'testing/output/baseline_methods_FDR_list.rds')

effect_info <- effect_info %>% right_join(rel_length, by = 'bandwidth') %>%
  mutate(across(`relative length`, ~ round(., 3))) %>% 
  dplyr::rename(`signal strength` = amplitude)

for(method in names(diff.testing.methods)){
  message(method)
  
  combined_fdr <- as_tibble(fdr.list[[method]]) %>% 
    mutate(gene = genes) %>%
    left_join(effect_info, by = 'gene')
  
  combined_fdr$jaccard <- sapply(seq_len(nrow(combined_fdr)), function(i) {
    compute_jaccard_sim(
      fdr = unlist(combined_fdr[i, clusters], use.names = FALSE),
      cell_type_range = cell_type_range,
      true_l = combined_fdr$interval_left[i],
      true_r = combined_fdr$interval_right[i],
      fdr_threshold = .05
    )
  })
  p <- combined_fdr %>% dplyr::rename(RL = `relative length`) %>%
        ggplot(aes(x = factor(`signal strength`), 
                               y = jaccard, 
                               fill = factor(`signal strength`))) +
    geom_boxplot(outlier.size = .75, alpha = 0.65) +
    facet_wrap(~ RL, labeller = label_both) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x    = method,
      y    = 'jaccard index', 
      fill = "signal strength"
    ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
      text            = element_text(family = "Arial", size = 6),
      axis.title      = element_text(family = "Arial", size = 8),
      axis.text       = element_text(family = "Arial", size = 6),
      strip.text      = element_text(family = "Arial", size = 8),
      legend.title    = element_text(family = "Arial", size = 8),
      legend.text     = element_text(family = "Arial", size = 6),
      panel.grid      = element_blank(),
      legend.position = "bottom",
      legend.box.spacing = unit(2, "pt"), 
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )
  
  print(p)
  ggsave(plot = p, filename = glue('testing/figs/{method}_JACCARD_20.pdf'), 
         width = 3, height = 2, units = 'in')    

}

