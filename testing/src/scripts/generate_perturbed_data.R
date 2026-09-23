generate_perturbed_data <- function(sce_fits, iter, G = 100, 
                                    perturb_frac = 0.2, N_SAMPLES = 3, 
                                    n_cores = 20, dataset= 'PANCREAS'){
  
  #' get learned parameter fits for the dataset from scDesign. Use baseline fitted mean parameters to introduce 
  #' localized perturbations along pseudotime across conditions. Use new parameters to simulate ground truth counts 
  #' data for differential testing. Save generated counts and metadata to file, return ground-truth parameter for perturbations to main
  #' 
  #' @param sce_fits list of fitted scDesign objects to use for simulations
  #' @param iter iteration number for unique storage of generated data
  #' @param G number of genes to include in perturbation experiment
  #' @param perturb_frac fraction of genes that must have perturbations
  #' @param N_SAMPLES number of randomly sampled replicates for each condition 
  #' @param n_cores number of cores used for simulating counts with scDesign3
  #' @return details of sampled genes and perturbations applied to each gene
  
  message(iter)
  set.seed(148 + iter)
  c(sce, sce.dat, sce.copula, sce.para) %<-% sce_fits #unlist learned variables for dataset
  
  ## sample genes to include in experiment
  genes <- sample(colnames(sce.dat$count_mat), size = G, replace = F); 
  t_ <- sce.dat$newCovariate$pseudotime
  
  ## extract the mean parameters for these genes
  mu_A <- sce.para$mean_mat[, genes]
  #====================================
  ## identify null and effect genes
  #====================================
  ord <- sample(G, replace = F)
  n_effect <- round(perturb_frac * G)
  effect_genes <- genes[ord[1:n_effect]]
  null_genes <- genes[ord[(n_effect + 1):G]]
  
  lfc_mat <- matrix(
    0,
    nrow = nrow(mu_A),
    ncol = G,
    dimnames = dimnames(mu_A)
  )
  
  #====================================
  # compute effects
  #====================================
  effect_info <- data.frame(
    gene = effect_genes,
    amplitude = sample(c(0.1, .25, .5, .75, 1, 1.5), size = length(effect_genes), replace = T), # sample discretely
    sign = sample(c(-1, 1), size = length(effect_genes), replace = T),
    center = runif(length(effect_genes), 0.2, 0.8),
    bandwidth = sample(
      c(0.025, 0.05, 0.075, 0.10, 0.15),
      length(effect_genes),
      replace = TRUE
    )
  )
  
  ## inject perturbations on mean parameter
  for (j in seq_len(effect_genes)) {
    g <- effect_info$gene[j]
    lfc_mat[, g] <-
      effect_info$sign[j] * effect_info$amplitude[j] *
      exp(
        -(t_ - effect_info$center[j])^2 /
          (2 * effect_info$bandwidth[j]^2)
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
             copula_list = sce.copula$copula_list$`1`[genes, genes],
             n_cores = n_cores,
             family_use = "nb",
             input_data = sce.dat$dat, # this is feature vector
             new_covariate = sce.dat$newCovariate,
             filtered_gene = sce.dat$filtered_gene)
  }
  
  counts_A <- vector("list", N_SAMPLES)
  counts_B <- vector("list", N_SAMPLES)
  
  for(r in seq_len(N_SAMPLES)){
    set.seed(148 + iter + r)
    counts_A[[r]] <- generate_counts(mu_A)
    counts_B[[r]] <- generate_counts(mu_B)
  }
  
  
  #===================================================================
  ## save counts matrix to memory
  #===================================================================
  sim_counts <- do.call(cbind, c(counts_A, counts_B))
  ## append unique sample name for each datapoint
  sim_meta <- data.frame(condition = factor(rep(c('condition_A', 'condition_B'), 
                                                each = N_SAMPLES * ncol(sce))),
                         replicate = rep(factor(rep(paste('replicate', 1:N_SAMPLES, sep = '_'), 
                                                    each = ncol(sce))), 
                                         2), 
                         sample = rep(factor(paste('sample', LETTERS[1:6], sep = '-')), 
                                      each = ncol(sce)),
                         pseudotime = rep(sce.dat$dat$pseudotime, 2*N_SAMPLES)) %>% as_tibble()
  
  
  write_rds(list(counts = sim_counts, meta = sim_meta),
            file = 'testing/output/{dataset}_iter_{iter}.rds')
  
  #===================================================================
  ## return list of experiment parameters 
  #===================================================================
  return(list(
        run_iter = iter,
        sampled_genes = genes, 
        null_genes = null_genes,
        effect_info = effect_info
        )
    )
  
}