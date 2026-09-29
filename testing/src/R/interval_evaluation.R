## evaluate argo

# run_interval_evaluation <- function(N = 50){
  message(glue('Running interval testing for Argo'))
  # read simulated data for this experiment iteration and unload the list into variables
  zeallot::`%<-%`(c(sim_counts, sim_meta, genes, effect_info), read_rds(glue('testing/data/temp/{dataset}_interval_perturbation.rds')))
  effect_info$gene <- make.unique(effect_info$gene)
  
  # get metadata for different scales of data representation
  zeallot::`%<-%`(c(cell_level_meta, bin_level_meta, sample_level_meta, .), prepare_metadata(sim_meta, N))
  
  ## convert to pseudobulk counts
  binned_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_bin)
  binned_counts <- varianceStabilizingTransformation(t(binned_counts))
  argo.input <- bin_level_meta %>%
    bind_cols(scale(t(binned_counts), center = T, scale = F))
  
  # map effect intervals into bin ids at experiment resolution for benchmarking
  coord_x <- unique(bin_level_meta$coord_x)
  ## append index of closest bins contained within the interval (can easily translate to containing the interval by moving by one unit)
  effect_info <- effect_info %>% mutate(
    lower_bin = findInterval(effect_info$interval_left, c(0, coord_x, 1)), # shift this back by 1 to get spanning bin
    upper_bin = findInterval(effect_info$interval_right, c(0, coord_x, 1)) - 1,  # don't include shift to get spanning bin
    bin_width = upper_bin - lower_bin + 1
  ) 
  
  boxplot(effect_info$bin_width~ effect_info$bandwidth)
  
  env <- environment() 
  
  zeallot::`%<-%`(c(res.argo, interval.stats, fits), 
                  run_argo(env, bin_meta = argo.input, 
                           FIT_CONDITIONAL_NOISE = TRUE, 
                           TEST_INTERVALS = TRUE))
  
  write_rds(fits, 'testing/output/argo_conditional_noise_interval_eval_fits.rds')
  
  pheatmap(interval.stats, cluster_cols = FALSE, show_rownames = FALSE, use_raster = FALSE,
           name = 'deviationn intervals')
  
  
  
  eval.metrics <- vector("list", length = nrow(effect_info))
  for(i in seq_len(nrow(effect_info))){
    eval.metrics[[i]] <- compute_overlap_statistics(intA = which(interval.stats[i, ] != 0), 
                                       intB = seq(effect_info$lower_bin[i], effect_info$upper_bin[i]))
  }
  
  effect_info <- right_join(effect_info, bind_rows(eval.metrics) %>% mutate(gene = effect_info$gene))
  # effect_info[names(eval.metrics[[1]])] <- bind_rows(eval.metrics)
  
  for(metric in names(eval.metrics[[1]])){
    message(metric)
    p <- ggplot(effect_info, aes(x = factor(amplitude), y = !!sym(metric), fill = factor(amplitude))) +
      geom_boxplot(outlier.size = 1.5, alpha = 0.8) +
      facet_wrap(~ bandwidth, labeller = label_both) + # Splits plots by bandwidth levels
      scale_fill_brewer(palette = "Set2") + # Professional color palette
      labs(
        title = glue("{metric} across Amplitude and Bandwidth"),
        x     = "Amplitude",
        y     = metric, 
        fill  = "Amplitude"
      ) +
      theme_bw(base_size = 16) + # Structured border theme for faceted plots
      theme(
        panel.grid      = element_blank(),
        legend.position    = "bottom"
      )
    
    ggsave(plot = p, filename = glue('testing/figs/{metric}_statistics_conditional_noise.pdf'), width = 8, height = 6)
  }
  
  
  
# }

