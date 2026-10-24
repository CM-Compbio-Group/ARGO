generate_perturbed_data <- function(expt = 1, expt_df = tibble(), sce_fits = list(), 
                                    N_GENES = 1e3, perturb_frac = 0.2, N_SAMPLES = 1, 
                                    dataset = 'PANCREAS'){
  
  #' get learned parameter fits for the dataset from scDesign. Use baseline fitted mean parameters to introduce 
  #' localized perturbations along pseudotime across conditions. Use new parameters to simulate ground truth counts 
  #' data for differential testing. 
  #' Save generated counts and metadata to file, return ground-truth perturbation parameters to main
  #' 
  #' @param expt experiment number for saving data
  #' @param expt_df perturbation plan
  #' @param sce_fits fitted scDesign3 parameters for the dataset
  #' @param N_GENES number of genes to include from data
  #' @param perturb_frac what fraction of the genes must have effects
  #' @param N_SAMPLES number of replicates to include for each condition
  #' @param dataset name of the dataset being studied
  #' 
  #' saves simulated data to rds file 
  #' script is hard coded for simulating only 2 conditions
  
  
  if(!(is.integer(expt)) && (!(is.numeric(expt)) ||  !(expt %% 1 == 0))){
    stop('Type Error : expt must be integer')
  }
  
  if(expt < 1){
    stop('Value Error : expt should be positive')
  }
  
  if(expt > nrow(expt_df)){
    stop('Value Error : expt should identify a row of expt_df')
  }
  
  stopifnot(
    "expt_df is missing or not a data frame" = is.data.frame(expt_df), 
    "expt_df does not have the required information" = all(c('amplitude', 'bandwidth') %in% (colnames(expt_df)))
  )
  
  stopifnot("input does not have the right objects" = identical(names(sce_fits), c('sce', 'data', 'copula', 'para')))
  
  if(!(is.double(N_GENES)) ||  !(N_GENES %% 1 == 0)){
    stop('Type Error : N_GENES must be integer')
  }
  if(N_GENES < 1){
    stop('Value Error : data must include at least one gene')
  }
  if(N_GENES > 2e4){
    stop('Value Error : N_GENES seems inflated. Check input.')
  }
  if(N_GENES < 10){
    warning('Consider increasing N_GENES to evaluate perturbation effects')
  }
  
  
  if(!(is.numeric(perturb_frac))){
    stop('Type Error : N_GENES must be numeric')
  }
  if(perturb_frac < 0 || perturb_frac > 1){
    stop('Value Error : perturb_frac must be between 0 and 1')
  }
  
  if(!(is.double(N_SAMPLES)) ||  !(N_SAMPLES %% 1 == 0)){
    stop('Type Error : N_SAMPLES must be integer')
  }
  if(N_SAMPLES < 1){
    stop('Value Error : N_SAMPLES must be positive')
  }
  if(N_SAMPLES > 100){
    stop('Value Error : N_SAMPLES seems inflated. Check input.')
  }
  
  # have to include bootstrapping for desired number of genes
  # include reasonable defaults
  # check columns of expt_df
  
  message(glue('simulating counts for experiment {expt}'))
  # set.seed(148 + expt)
  
  amplitude = expt_df$amplitude[expt]
  bandwidth = expt_df$bandwidth[expt]
  
  
  zeallot::`%<-%`(c(sce, sce.dat, sce.copula, sce.para), sce_fits) #unlist learned variables for dataset
  
  ## sample genes to include in experiment
  genes <- sample(colnames(sce.dat$count_mat), size = N_GENES, replace = T); 
  unique_genes <- make.unique(genes) # make gene names unique when we have repeated samples
  t_ <- sce.dat$newCovariate$pseudotime
  N_TOTAL <- 2*N_SAMPLES # total number of samples simulated
  
  ## get fitted parameters for the sampled genes
  mu_A <- sce.para$mean_mat[, genes]
  sigma <- sce.para$sigma_mat[, genes]
  zero <- sce.para$zero_mat[, genes]
  colnames(mu_A) <- colnames(sigma) <- colnames(zero) <- unique_genes
  
  copula <- sce.copula$copula_list$`1`
  copula <- copula[genes, genes] # get gene-gene correlations for bootstrapped genes
  dimnames(copula) <- list(unique_genes, unique_genes)
  copula_list <- list(`1` = copula)
  
  ## update single-cell experiment object
  sce <- sce[genes, ]
  rownames(sce) <- unique_genes
  #====================================
  ## identify null and effect genes
  #====================================
  ord <- sample(N_GENES, replace = F)
  n_effect <- round(perturb_frac * N_GENES)
  effect_genes <- unique_genes[ord[1:n_effect]]
  null_genes <- unique_genes[ord[(n_effect + 1):N_GENES]]
  
  lfc_mat <- matrix(
    0,
    nrow = nrow(mu_A),
    ncol = N_GENES,
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
    simu_new(sce = sce,
             mean_mat = mu,
             sigma_mat = sigma,
             zero_mat = zero,
             quantile_mat = NULL,
             copula_list = copula_list,
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
            file = glue('testing/data/10_20/{dataset}_expt_{expt}.rds'))
  
}
