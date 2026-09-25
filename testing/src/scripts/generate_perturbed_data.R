generate_perturbed_data <- function(amplitude, bandwidth, expt, env){
  
  #' get learned parameter fits for the dataset from scDesign. Use baseline fitted mean parameters to introduce 
  #' localized perturbations along pseudotime across conditions. Use new parameters to simulate ground truth counts 
  #' data for differential testing. Save generated counts and metadata to file, return ground-truth parameter for perturbations to main
  #' 
  #' @param amplitude magnitude of perturbation FC 
  #' @bandwidth extent of perturbation between conditions
  #' @param expt experiment number for saving data
  #' @param env parent environment for scoping callstack
  #' 
  #' saves simulated data to rds file 
  #' script is hard coded for simulating only 2 conditions
  
  message(glue('simulating counts for experiment {expt}'))
  set.seed(148 + expt)
  c(sce, sce.dat, sce.copula, sce.para) %<-% env$sce_fits #unlist learned variables for dataset
  
  ## sample genes to include in experiment
  genes <- sample(colnames(sce.dat$count_mat), size = env$G, replace = F); 
  t_ <- sce.dat$newCovariate$pseudotime
  N_TOTAL <- 2*env$N_SAMPLES # total number of samples simulated
  
  ## extract the mean parameters for these genes
  mu_A <- sce.para$mean_mat[, genes]
  sce.copula$copula_list$`1` <- sce.copula$copula_list$`1`[genes, genes]
  #====================================
  ## identify null and effect genes
  #====================================
  ord <- sample(env$G, replace = F)
  n_effect <- round(env$perturb_frac * env$G)
  effect_genes <- genes[ord[1:n_effect]]
  null_genes <- genes[ord[(n_effect + 1):env$G]]
  
  lfc_mat <- matrix(
    0,
    nrow = nrow(mu_A),
    ncol = env$G,
    dimnames = dimnames(mu_A)
  )
  #====================================
  # compute effects
  #====================================
  # amplitude and bandwidth are fixed, only vary where the perturbation occurs, and sign 
  effect_info <- data.frame(
    gene = effect_genes,
    amplitude = amplitude, 
    bandwidth = bandwidth,
    sign = sample(c(-1, 1), size = length(effect_genes), replace = T),
    center = runif(length(effect_genes), 0.2, 0.8)
    )
  
  
  ## inject perturbations on mean parameter
  for (j in seq_along(effect_genes)) {
    g <- effect_info$gene[j]
    lfc_mat[, g] <-
      effect_info$sign[j] * amplitude *
      exp(
        -(t_ - effect_info$center[j])^2 /
          (2 * bandwidth^2)
      )
  }
  
  ## compute the new mean parameter with these effects
  mu_B <- mu_A * exp(lfc_mat)
  #===================================================================
  ## simulate counts data
  #===================================================================
  generate_counts <- function(mu){
    simu_new(sce = sce[genes, ],
             mean_mat = mu,
             sigma_mat = sce.para$sigma_mat[, genes],
             zero_mat = sce.para$zero_mat[, genes],
             quantile_mat = NULL,
             copula_list = sce.copula$copula_list,
             n_cores = env$n_cores,
             family_use = "nb",
             input_data = sce.dat$dat, # this is feature vector
             new_covariate = sce.dat$newCovariate,
             filtered_gene = sce.dat$filtered_gene)
  }
  
  counts_A <- vector("list", env$N_SAMPLES)
  counts_B <- vector("list", env$N_SAMPLES)
  
  for(r in seq_len(env$N_SAMPLES)){
    counts_A[[r]] <- generate_counts(mu_A)
    counts_B[[r]] <- generate_counts(mu_B)
  }
  
  
  #===================================================================
  ## save counts matrix to memory
  #===================================================================
  sim_counts <- do.call(cbind, c(counts_A, counts_B))
  ## append unique sample name for each datapoint
  sim_meta <- data.frame(condition = factor(rep(c('condition_A', 'condition_B'), 
                                                each = env$N_SAMPLES * ncol(sce))),
                         replicate = rep(factor(rep(paste('replicate', 1:env$N_SAMPLES, sep = '_'), 
                                                    each = ncol(sce))), 
                                         2), 
                         sampleID = rep(factor(paste('sample', LETTERS[1:N_TOTAL], sep = '-')), 
                                      each = ncol(sce)),
                         pseudotime = rep(sce.dat$dat$pseudotime, N_TOTAL)) %>% as_tibble()
  
  
  write_rds(list(counts = sim_counts, 
                 meta = sim_meta, 
                 effect_genes = effect_genes, 
                 null_genes = null_genes),
            file = glue('testing/data/temp/{env$dataset}_expt_{expt}.rds'))
  
}

#===================================================================
## return list of experiment parameters 
#===================================================================
# return(list(
#       run_iter = expt,
#       sampled_genes = genes, 
#       null_genes = null_genes,
#       effect_info = effect_info
#       )
#   )