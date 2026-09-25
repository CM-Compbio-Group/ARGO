#' use scDesign3 to simulate ground truth single cell expression data along trajectories
#' that can be used to benchmark argo against baseline methods for differential testing.
#' We start by using scDesign to fit parameters for each gene in the dataset
#' We then perturb the mean parameters along the trajectory to create different kinds of effects, 
#' and check the capacity of algorithms to recover these deviances

rm(list=ls())
source('testing/src/scripts/00_load_packages_and_functions.R')
dataset <- 'PANCREAS'

## make a map of parameter values and experiments
perturb_amplitude <- c(0.3, 0.6, 1) # magnitude of LFC in perturbation
perturb_bandwidth <- c(0.03, 0.06, 0.1) # width of the perturbation
N_EXPERIMENTS <- 5 # no. of replicates for each parameter combination
expt_df <- crossing(amplitude = perturb_amplitude,
                    bandwidth = perturb_bandwidth,
                    iter = factor(paste('iter', seq(N_EXPERIMENTS), sep = '-'))) %>% 
  mutate(Experiment = row_number())



FIT_DATA <- FALSE
GENERATE_DATA <- FALSE

if(FIT_DATA){
  sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))
  print(sce) ## 1000 genes, 2087 cells
  sce_fits <- fit_scDesign3(sce) # learn mean and variance parameters for all genes
  write_rds(sce_fits, glue('testing/output/{dataset}_sce_fits.rds'))
}




if(GENERATE_DATA){
  ## define control variables for perturbation experiment
  G = 500 # number of genes to use
  perturb_frac = 0.2 # fraction of genes that must have perturbations
  N_SAMPLES = 3 # number of randomly sampled replicates for each condition 
  n_cores = 20 # number of cores used for simulating counts with scDesign3
  env <- environment()
  
  sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))
  
  message('generating perturbation datasets')
  for(j in seq_len(nrow(expt_df))){
    message(glue('Expt-{j}'))
    generate_perturbed_data(amplitude = expt_df$amplitude[j], 
                            bandwidth = expt_df$bandwidth[j], 
                            expt = j, env = env)
  }
}


## test on one experiment
# metrics <- run_differential_testing(expt = 20, dataset)
# summary(metrics)

## run differential testing and gather evaluation metrics on each dataset
message('now running differential testing')
start_time <- Sys.time()
eval.metrics <- lapply(seq_len(nrow(expt_df)), 
  run_differential_testing, dataset, N = 50, M = 100)
end_time <- Sys.time()
print(end_time - start_time)

write_rds(eval.metrics, 'testing/output/diff_testing_results_EXPT3.rds')


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

# dplyr::filter(eval_long, method == 'argo_conditional_noise' & amplitude == 0.3 & bandwidth == 0.03 & metric == 'F1')
# View(filter(eval_long, Experiment == 10))


for(metrx in c('AUPRC', 'F1')){
  message(metrx)
  p <- ggplot(eval_long %>% filter(metric == metrx), aes(x = method, y = value, color = method)) +
    # Scales the 10 iterations tightly within each of the 7 method columns
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
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
  
  
  ggsave(plot = p, filename = glue('testing/figs/{metrx}_comparative_EXPT2.pdf'), 
         device = 'pdf', width = 15, height = 15)

}

