generate_interval_perturbation <- function(sce_fits, dataset){
  
  #' perturbing all genes with different combinations of amplitude and bandwidth parameters
  #' 
  #' @param sce_fits fitted scDesign3 parameters for the dataset
  #' @param dataset name of the dataset being studied
  #' 
  #' saves simulated data to rds file 
  #' script is hard coded for simulating only 2 conditions
  
  set.seed(148)
  zeallot::`%<-%`(c(sce, sce.dat, sce.copula, sce.para), sce_fits) #unlist learned variables for dataset
  
  N_SAMPLES <- 1
  N_TOTAL <- 2*N_SAMPLES
  ## sample genes to include in experiment
  effect_genes <- sample(rownames(sce), size = 3000, replace = T)
  unique_effect_genes <- make.unique(effect_genes)
  
  t_ <- sce.dat$newCovariate$pseudotime
  ord <- order(t_)
  
  ## extract the mean parameters for these genes
  mu_A <- sce.para$mean_mat[, effect_genes]
  sigma_A <- sce.para$sigma_mat[, effect_genes]
  zero_A <- sce.para$zero_mat[, effect_genes]
  
  colnames(mu_A) <- colnames(sigma_A) <- colnames(zero_A) <- unique_effect_genes
  ## update the correlation matrix b/w genes
  # sce.copula$copula_list$`1` <- sce.copula$copula_list$`1`[effect_genes, effect_genes]
  copula <- sce.copula$copula_list$`1`
  copula <- copula[effect_genes, effect_genes]
  dimnames(copula) <- list(unique_effect_genes, unique_effect_genes)
  copula_list <- list(`1` = copula)
  sce <- sce[effect_genes, ]
  rownames(sce) <- unique_effect_genes
  
  #====================================
  ## identify null and effect genes
  #====================================
  lfc_mat <- matrix(
    0,
    nrow = nrow(mu_A),
    ncol = ncol(mu_A),
    dimnames = dimnames(mu_A)
  )
  
  #====================================
  # compute effects
  #====================================
  # amplitude and bandwidth are fixed, only vary where the perturbation occurs, and sign 
  effect_info <- tibble(
    gene = unique_effect_genes,
    amplitude = sample(c(0.3, 0.6, 1), size = length(effect_genes), replace = T), 
    bandwidth = sample(c(0.03, 0.06, 0.1), size = length(effect_genes), replace = T),
    sign = sample(c(-1, 1), size = length(effect_genes), replace = T),
    center = runif(length(effect_genes), 0.2, 0.8)
  )
  
  print(table(effect_info$amplitude, effect_info$bandwidth))
  
  ## inject perturbations on mean parameter
  for (j in seq_along(effect_genes)) {
    lfc_mat[, j] <-
      effect_info$sign[j] * effect_info$amplitude[j] *
      exp(
        -(t_ - effect_info$center[j])^2 /
          (2 * effect_info$bandwidth[j]^2)
      )
  }
  
  ## visualize the perturbation map
  # pheatmap(lfc_mat[ord, ], use_raster=FALSE, cluster_rows = F, show_rownames = F, show_colnames = F, name = 'effect size')
  
  ## for every gene, store the interval with non-zero deviation
  thresh = 0.1
  lower_time_bound <- apply(lfc_mat, 2, function(z) min(t_[abs(z) > thresh])) # min time point where absolute lfc exceeds threshold
  upper_time_bound <- apply(lfc_mat, 2, function(z) max(t_[abs(z) > thresh]) ) # upper bound of interval
  effect_info <- effect_info %>% mutate(interval_left_01 = lower_time_bound, 
                         interval_right_01 = upper_time_bound)
  
  ## compute the new mean parameter with these effects
  mu_B <- mu_A * exp(lfc_mat)
  #===================================================================
  ## simulate counts data
  #===================================================================
  generate_counts <- function(mu){
    simu_new(sce = sce,
             mean_mat = mu,
             sigma_mat = sigma_A,
             zero_mat = zero_A,
             quantile_mat = NULL,
             copula_list = copula_list,
             n_cores = 20,
             family_use = "nb",
             input_data = sce.dat$dat, # this is feature vector
             new_covariate = sce.dat$newCovariate,
             filtered_gene = sce.dat$filtered_gene)
  }
  
  counts_A <- generate_counts(mu_A)
  counts_B <- generate_counts(mu_B)
  
  #===================================================================
  ## save counts matrix to memory
  #===================================================================
  ## append unique sample name for each cell
  sim_meta <- tibble(cell_id = rep(colnames(sce), N_TOTAL),
                         condition = factor(rep(c('condition_A', 'condition_B'), 
                                                each = N_SAMPLES * ncol(sce))),
                         sampleID = rep(factor(paste('sample', LETTERS[1:N_TOTAL], sep = '-')), 
                                        each = ncol(sce)),
                         pseudotime = rep(sce.dat$dat$pseudotime, N_TOTAL), 
                         cell_type = rep(sce$cell_type, N_TOTAL)) %>% as_tibble()
  
  sim_meta <- sim_meta %>% mutate(cell_id = paste(sampleID, cell_id, sep = '_'))
  sim_counts <- as(cbind(counts_A, counts_B), 'sparseMatrix')
  colnames(sim_counts) <- sim_meta$cell_id
  stopifnot(!any(is.na(sim_counts)))
  print(glue("size of generated counts data is {nrow(sim_counts)}"))
  
  write_rds(list(counts = sim_counts, 
                 meta = sim_meta, 
                 genes = effect_genes, 
                 effect_info = effect_info),
            file = glue('testing/data/temp/{dataset}_interval_perturbation.rds'))
  
}
