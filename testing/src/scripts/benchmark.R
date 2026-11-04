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
  N_SAMPLES = 3 # number of replicates per condition 
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
                               run_differential_testing, N = 25, dataset) |> futurize(seed = TRUE)
  end_time <- Sys.time()
  print(end_time - start_time)
  write_rds(eval.metrics, glue('testing/output/diff_testing_N_25_p_500_eval_metrics.rds'))
  plan(sequential)
}

eval <- run_differential_testing(expt=29, N=25)

## find the relative length that the 20% threshold corresponds to



if(GENERATE_INTERVAL_BENCHMARK){
  set.seed(148)
  sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))
  generate_interval_perturbation(sce_fits = sce_fits, N_GENES = 3e3, N_SAMPLES = 3, dataset)
}

if(RUN_INTERVAL_BENCHMARK){
  run_interval_evaluation( N = 25, dataset = dataset)  
}
effect_info <- read_rds('testing/output/interval_effect_info.rds')

## find mean eval metric by parameter combination
effect_info %>% group_by(amplitude, bandwidth) %>% summarise(mean(jaccard_index, na.rm = T))



if(PLOT_RESULTS){
  eval.metrics <- read_rds(glue('testing/output/diff_testing_N_25_p_500_eval_metrics.rds'))
  metrics <- do.call(rbind, eval.metrics)
  method <- rownames(metrics)
  
  eval_long <- metrics %>% as_tibble() %>% 
    mutate(method = factor(method, levels = c('DESeq2', 'edgeR', 'Limma', 'tradeSeq', 
                                              'Lamian', 'argo_conditional_noise')), 
           Experiment = rep(seq(nrow(expt_df)), each = 6)) %>%
    pivot_longer(
      cols = c(sensitivity, specificity, precision, empirical_FDR, accuracy, F1, AUPRC),
      names_to = "metric",
      values_to = "value"
    ) %>% merge(expt_df, by = 'Experiment') %>% 
    mutate(method = fct_recode(method, 
                               'ARGO' = 'argo_conditional_noise')) 
  
  
  eval_long <- left_join(eval_long, rel_length, by = 'bandwidth') %>% 
               mutate(across(rel_length, ~ round(., 3))) %>% 
               dplyr::rename(`relative length` = rel_length,
                                    `signal strength` = amplitude)
  
  method_cols <- c(
    "ARGO"     = "orchid4",
    "DESeq2"   = "#bebada",
    "edgeR"    = "#8dd3c7",
    "Limma"    = "lightskyblue2",
    "Lamian"   = "lightgreen",
    "tradeSeq" = "lightpink2"
  )
  
  for(metrx in c('AUPRC', 'F1')){
    message(metrx)
    p <- ggplot(eval_long %>% filter(metric == metrx) %>%
                  mutate(`relative length` = scales::percent(`relative length`)),
                aes(x = method, y = value, color = method)) +
      geom_jitter(width = 0.125, alpha = 0.7, size = 1.25) + 
      stat_summary(fun = mean, geom = "point", color = "black", size = 2, pch = 18) +
      facet_wrap(
        ~ `signal strength` * `relative length`, 
        scales = "free_y", ncol = 3,
        labeller = \(x) label_both(x, multi_line = TRUE, sep = ": ")
      ) +
      theme_minimal(base_family = "Arial", base_size = 6) +
      scale_color_manual(values = method_cols) +
      scale_x_discrete(
        # 1. Shrink outer padding around the first and last x categories
        expand = expansion(add = c(0.3, 0.3)),
        labels = function(x) str_wrap(x, width = 10)
      ) +
      labs(x = NULL, y = metrx) +
      theme(
        panel.grid = element_blank(),
        
        # Add border around each subpanel to create a grid between them
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
        
        # 2. Controls panel width-to-height ratio: taller panels pull x-coordinates closer together
        aspect.ratio = 1.1,
        
        # 3. Tighten distance from plot border to x tick labels
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6, 
                                   margin = margin(t = 2, unit = "pt")),
        axis.ticks.length.x = unit(2, "pt"),
        
        # 4. Tighten y-axis text/tick margins to minimize the gap between facet columns
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 6, 
                                   margin = margin(r = 1, unit = "pt")),
        axis.ticks.length.y = unit(2, "pt"),
        
        legend.position = "none",
        strip.text = element_text(face = "plain", size = 6, margin = margin(t = 1, b = 1)),
        
        # 5. Inter-panel spacing
        panel.spacing.x = unit(5, "pt"),
        panel.spacing.y = unit(5, "pt"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    )
    
    print(p)
    # Reduce output width from 5 to 4 to further compress x spacing across the 3 columns
    ggsave(
      plot = p, 
      filename = glue('testing/figs/{metrx}_evaluation.pdf'), 
      device = 'pdf', 
      width = 4, 
      height = 4.5, 
      units = 'in'
    )
    
  }
  
}