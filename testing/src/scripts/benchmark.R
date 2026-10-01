#' use scDesign3 to simulate ground truth single cell expression data along trajectories
#' that can be used to benchmark argo against baseline methods for differential testing.
#' We start by using scDesign to fit parameters for each gene in the dataset
#' We then perturb the mean parameters along the trajectory to create different kinds of effects, 
#' and check the capacity of algorithms to recover these deviances
#' 
#' Adding benchmarks for the deviation interval recognition property of ARGO -- simulating genes with different 
#' perturbation amplitudes and bandwidths, and comparing recovered bands with ground truth regions with non-zero FC


rm(list=ls())
source('testing/src/scripts/00_load_packages_and_functions.R')
dataset <- 'PANCREAS'


FIT_DATA <- TRUE
GENERATE_DATA <- TRUE
RUN_TESTING <- TRUE
PLOT_RESULTS <- TRUE

GENERATE_INTERVAL_BENCHMARK <- TRUE
RUN_INTERVAL_BENCHMARK <- TRUE


## make a map of parameter values and experiments
N_EXPERIMENTS <- 5 # no. of replicates for each parameter combination
perturb_amplitude <- c(0.3, 0.6, 1) # magnitude of LFC in perturbation
perturb_bandwidth <- c(0.03, 0.06, 0.1) # width of the perturbation
expt_df <- crossing(amplitude = perturb_amplitude,
                    bandwidth = perturb_bandwidth,
                    iter = factor(paste('iter', seq(N_EXPERIMENTS), sep = '-'))) %>% 
  mutate(Experiment = row_number())


if(FIT_DATA){
  sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))
  print(sce) ## 1000 genes, 2087 cells
  sce_fits <- fit_scDesign3(sce) # learn mean and variance parameters for all genes
  write_rds(sce_fits, glue('testing/output/{dataset}_sce_fits.rds'))
}



if(GENERATE_DATA){
  ## define control variables for perturbation experiment
  N_SAMPLES = 3 # number of randomly sampled replicates for each condition 
  N_GENES = 500 # number of genes to use
  perturb_frac = 0.2 # fraction of genes that must have perturbations
  sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))
  
  ## generate the datasets in parallel 
  # generate_perturbed_data(expt=29, expt_df, sce_fits,
  #           N_GENES = 2000, perturb_frac = 0.2, N_SAMPLES = 1, dataset) # test run 
  message('generating perturbation datasets')
  base::lapply(seq_len(nrow(expt_df)), generate_perturbed_data, expt_df, sce_fits,
         N_GENES, perturb_frac, N_SAMPLES, dataset) |> futurize(seed = TRUE)
}



## run differential testing and gather evaluation metrics on each dataset
if(RUN_TESTING){
  message('now running differential testing')
  plan(multisession, workers = 8)
  start_time <- Sys.time()
  eval.metrics <- base::lapply(seq_len(nrow(expt_df)), 
                               run_differential_testing, N = 50, dataset) |> futurize(seed = TRUE)
  end_time <- Sys.time()
  print(end_time - start_time)
  write_rds(eval.metrics, glue('testing/output/diff_testing_eval_metrics.rds'))
  plan(sequential)
}


if(GENERATE_INTERVAL_BENCHMARK){
  set.seed(148)
  sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))
  generate_interval_perturbation(sce_fits = sce_fits, N_GENES = 3e2, N_SAMPLES = 1, dataset)
}

if(RUN_INTERVAL_BENCHMARK){
  run_interval_evaluation(dataset, N = 50)  
}


if(PLOT_RESULTS){
  metrics <- do.call(rbind, eval.metrics)
  method <- rownames(metrics)
  eval_long <- metrics %>% as_tibble() %>% 
    mutate(method = method, 
           Experiment = rep(seq(nrow(expt_df)), each = 7)) %>%
    mutate(method = factor(method, levels = c('DESeq2', 'edgeR', 'Limma', 'tradeSeq', 'Lamian', 'argo_shared_noise', 'argo_conditional_noise'))) %>%
    pivot_longer(
      cols = c(sensitivity, specificity, precision, empirical_FDR, accuracy, F1, AUPRC),
      names_to = "metric",
      values_to = "value"
    )
  
  eval_long <- merge(eval_long, expt_df, by = 'Experiment') %>% as_tibble()
  eval_long
  
  for(metrx in c('AUPRC', 'F1')){
    message(metrx)
    p <- ggplot(eval_long %>% filter(metric == metrx), aes(x = method, y = value, color = method)) +
      # Scales the 10 iterations tightly within each of the 7 method columns
      geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) + 
      # ylim(c(0, 1)) + 
      # Adds a clear center-point for the mean
      stat_summary(fun = mean, geom = "point", color = "black", size=3, pch=18) +
      # Splitting into metric subpanels
      facet_wrap(~ amplitude * bandwidth, scales = "free_y", ncol = 3, 
                 labeller = \(x) label_both(x, multi_line = TRUE, sep = ": ")) +
      theme_minimal() +
      scale_color_brewer(palette = "Set2") + 
      labs(
        x = "Method",
        y = metrx,
        title = glue("Metric Performance Across {N_EXPERIMENTS} Iterations")
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(vjust = 1, hjust = 1, size = 12),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        panel.spacing.x = unit(2.5, "lines"),  # Horizontal space between columns
        panel.spacing.y = unit(2.5, "lines"),   # Increased vertical space between rows
        plot.margin = margin(t = 20, r = 20, b = 30, l = 20, unit = "pt") 
      ) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
    
    
    ggsave(plot = p, filename = glue('testing/figs/{metrx}_evaluation.pdf'), 
           device = 'pdf', width = 15, height = 15)
    
  }
}