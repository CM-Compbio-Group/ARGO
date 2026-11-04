run_interval_evaluation <- function(N = 25, FIT_CONDITIONAL_NOISE = TRUE, SAVE_FITS = TRUE, dataset = 'PANCREAS'){
  #' compare intervals of significant deviation evaluated by ARGO and compare with ground truth
  #' genes are perturbed with different strength (amplitude) and extent (bandwidth) of FC between conditions
  #' We use a range of metrics to evaluate how well the recovered intervals from ARGO recover domains of significant deviation
  #' @param N number of bins to split temporal trajectory into
  #' @param FIT_CONDITIONAL_NOISE whether to fit ARGO with condition-specific noise
  #' @param SAVE_FITS whether we want to save the ARGO differential testing results 
  #' @param dataset ground truth experiment perturbation data is coming from
  #' 
  #' @returns () : directly stores the ARGO fits and visualization of metrics in corresponding destinations
  #' 
  #' @error-checks check if input file exists, value and type errors for input parameters
  
  ## validate inputs
  if(!(is.double(N)) ||  !(N %% 1 == 0)){
    stop('Type Error : N must be integer')
  }
  if(!is.logical(FIT_CONDITIONAL_NOISE)){
    stop('Type Error : FIT_CONDITIONAL_NOISE must be boolean')
  }
  if(!is.logical(SAVE_FITS)){
    stop('Type Error : SAVE_FITS must be boolean')
  }
  if(N < 1 || N > 1000){
    stop('Value Error : Value of N is not suitable for pseudo-bulking')
  }
  if(N < 10){
    warning('Temporal resolution might not be sufficient for recovering perturbations. Consider increasing N')
  }
  if(N > 100){
    warning('Bins might be too sparse. Consider lowering N')
  }
  
  stopifnot("input file not found" = file.exists(glue('testing/data/temp/{dataset}_w_replicates_interval_perturbation.rds')))
  
  message(glue('Running interval recovery evaluation for Argo'))
  # read simulated data for this experiment iteration and unload the list into variables
  zeallot::`%<-%`(c(sim_counts, sim_meta, genes, effect_info), 
                  read_rds(glue('testing/data/temp/{dataset}_w_replicates_interval_perturbation.rds')))
  
  effect_info <- effect_info %>% 
    left_join(rel_length, by = 'bandwidth') %>% 
    dplyr::rename(`signal strength` = amplitude) %>%
    mutate(across(`relative length`, ~ round(., 3)))
  
  
  # get metadata for different scales of data representation
  zeallot::`%<-%`(c(cell_level_meta, bin_level_meta,
                    sample_level_meta, cell_type_meta), 
                  prepare_metadata(sim_meta, N))
  
  ## convert to pseudobulk counts
  binned_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_bin)
  binned_counts <- varianceStabilizingTransformation(t(binned_counts))
  argo.input <- bin_level_meta %>%
    bind_cols(scale(t(binned_counts), center = T, scale = T))
  # map effect intervals into bin ids at experiment resolution for benchmarking
  coord_x <- unique(bin_level_meta$coord_x)
  
  
  effect_info <- effect_info %>%
    mutate(
      lower_bin = findInterval(effect_info$interval_left, c(0, coord_x, 1)), # shift this back by 1 to get spanning bin
      upper_bin = findInterval(effect_info$interval_right, c(0, coord_x, 1)) - 1,  # don't include shift to get spanning bin
      bin_width = upper_bin - lower_bin + 1
    ) 
  ## append index of closest bins contained within the interval (can easily translate to containing the interval by moving by one unit)
  boxplot(effect_info$bin_width~ effect_info$bandwidth)
  
  env <- environment() 
  system.time({
    zeallot::`%<-%`(c(res.argo, interval.stats, fits), 
                    run_argo(env, bin_meta = argo.input, 
                             FIT_CONDITIONAL_NOISE = TRUE, 
                             TEST_INTERVALS = TRUE))
  })
  
  write_rds(fits, glue('testing/output/argo_fits_w_replicates_NBINS_{N}.rds'))
  
  # pheatmap(interval.stats, cluster_cols = FALSE, show_rownames = FALSE, use_raster = FALSE,
  #          name = 'deviationn intervals')
  
  # fits <- read_rds(glue('testing/output/argo_fits_w_replicates_NBINS_{25}.rds'))
  eval.metrics <- vector("list", length = nrow(effect_info))
  for(i in seq_along(genes)){
    eval.metrics[[i]] <- compute_overlap_statistics(intA = which(interval.stats[i, ] != 0), 
                                       intB = seq(effect_info$lower_bin[i], effect_info$upper_bin[i]))
  }
  
  eval.metrics <- bind_rows(eval.metrics) %>% 
    mutate(gene = effect_info$gene) %>%
    left_join(effect_info, by = 'gene')
  
  write_rds(eval.metrics, file = glue('testing/output/argo_interval_effect_info_NBINS_{N}.rds'))
  
  for(metric in colnames(eval.metrics[1:5])){
    message(metric)
    p <- ggplot(eval.metrics, aes(x = factor(`signal strength`), 
                                 y = !!sym(metric), 
                                 fill = factor(`signal strength`))) +
      geom_boxplot(outlier.size = 1.5, alpha = 0.8) +
      facet_wrap(~ `relative length`, labeller = label_both) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        x    = "ARGO",
        y    = metric, 
        fill = "signal strength"
      ) +
      theme_bw(base_size = 8, base_family = "Arial") +
      theme(
        text            = element_text(family = "Arial", size = 8),
        axis.title      = element_text(family = "Arial", size = 8),
        axis.text       = element_text(family = "Arial", size = 8),
        strip.text      = element_text(family = "Arial", size = 8),
        legend.title    = element_text(family = "Arial", size = 8),
        legend.text     = element_text(family = "Arial", size = 8),
        panel.grid      = element_blank(),
        legend.position = "bottom",
        legend.box.spacing = unit(2, "pt"), 
        
        # 2. Remove any top padding inside the legend container itself
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
      )
    
    ggsave(plot = p, filename = glue('testing/figs/{metric}_statistics_NBINS_{N}_signal_20.pdf'), 
           width = 4, height = 2.5, units = 'in')
    
  }
  
}

