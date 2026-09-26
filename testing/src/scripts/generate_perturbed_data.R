generate_perturbed_data <- function(expt, expt_df, sce_fits, G, perturb_frac, N_SAMPLES, dataset){
  
  #' get learned parameter fits for the dataset from scDesign. Use baseline fitted mean parameters to introduce 
  #' localized perturbations along pseudotime across conditions. Use new parameters to simulate ground truth counts 
  #' data for differential testing. Save generated counts and metadata to file, return ground-truth parameter for perturbations to main
  #' 
  #' @param expt experiment number for saving data
  #' @param expt_df perturbation plan
  #' @param sce_fits fitted scDesign3 parameters for the dataset
  #' @param G number of genes to include from data
  #' @param perturb_frac what fraction of the genes must have effects
  #' @param N_SAMPLES number of replicates to include for each condition
  #' @param dataset name of the dataset being studied
  #' 
  #' saves simulated data to rds file 
  #' script is hard coded for simulating only 2 conditions
  
  message(glue('simulating counts for experiment {expt}'))
  set.seed(148 + expt)
  
  amplitude = expt_df$amplitude[expt]
  bandwidth = expt_df$bandwidth[expt]
  
  
  zeallot::`%<-%`(c(sce, sce.dat, sce.copula, sce.para), sce_fits) #unlist learned variables for dataset
  
  ## sample genes to include in experiment
  genes <- sample(colnames(sce.dat$count_mat), size = G, replace = F); 
  t_ <- sce.dat$newCovariate$pseudotime
  N_TOTAL <- 2*N_SAMPLES # total number of samples simulated
  
  ## extract the mean parameters for these genes
  mu_A <- sce.para$mean_mat[, genes]
  sce.copula$copula_list$`1` <- sce.copula$copula_list$`1`[genes, genes]
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
             n_cores = 1,
             family_use = "nb",
             input_data = sce.dat$dat, # this is feature vector
             new_covariate = sce.dat$newCovariate,
             filtered_gene = sce.dat$filtered_gene)
  }
  
  counts_A <- vector("list", N_SAMPLES)
  counts_B <- vector("list", N_SAMPLES)
  
  for(r in seq_len(N_SAMPLES)){
    counts_A[[r]] <- generate_counts(mu_A)
    counts_B[[r]] <- generate_counts(mu_B)
  }
  
  
  #===================================================================
  ## save counts matrix to memory
  #===================================================================
  ## append unique sample name for each cell
  sim_meta <- data.frame(cell_id = rep(colnames(sce), N_TOTAL),
                         condition = factor(rep(c('condition_A', 'condition_B'), 
                                  each = N_SAMPLES * ncol(sce))),
                         replicate = rep(factor(rep(paste('replicate', 1:N_SAMPLES, sep = '_'), 
                                                    each = ncol(sce))), 2), 
                         sampleID = rep(factor(paste('sample', LETTERS[1:N_TOTAL], sep = '-')), 
                                      each = ncol(sce)),
                         pseudotime = rep(sce.dat$dat$pseudotime, N_TOTAL), 
                         cell_type = rep(sce$cell_type, N_TOTAL)) %>% as_tibble()
  
  
  sim_meta <- sim_meta %>% mutate(cell_id = paste(sampleID, cell_id, sep = '_'))
  sim_counts <- as(do.call(cbind, c(counts_A, counts_B)), 'sparseMatrix')
  colnames(sim_counts) <- sim_meta$cell_id
  
  write_rds(list(counts = sim_counts, 
                 meta = sim_meta, 
                 effect_genes = effect_genes, 
                 null_genes = null_genes, 
                 effect_info = effect_info),
            file = glue('testing/data/temp/{dataset}_expt_{expt}.rds'))
  
}
