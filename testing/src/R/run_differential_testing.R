run_differential_testing <- function(expt = 1, N = 25, dataset = 'PANCREAS'){
  #' load simulated counts for the iteration and run differential testing across all algorithms
  #' simulated counts are first converted into pseudo-bulk counts with desired number of bins along the isodepth
  #' raw / normalized pseudo-bulk counts are used as input for the different testing algorithms
  #' differential testing results are returned as a list to caller function
  #' 
  #' Run edgeR, DESeq2 and limma with cell type annotation-level pseudobulk counts data
  #' Run Lamian and tradeSeq with single-cell resolution data -- individual cell's expression and pseudotime values
  #' Run and test ARGO on 50 bins along pseudotime -- bins are chosen to have equal number of cells, so distance between them will vary
  #'  no need to compute the confidence bands with argo -- only need the global test statistic / p.value
  #' 
  #' @param expt experiment iteration
  #' @param dataset dataset label 
  #' @param N number of observation points along pseudotime
  #' Expression counts etc are read from memory, not passed as arguments: 
  #'  a list of the generated counts, metadata, and the set of effect and null genes
  #' 
  #' @returns list of evaluation metrics across differential testing methods on this simulation run
  #' 
  
  ## validate inputs
  if(!(is.double(expt)) ||  !(expt %% 1 == 0)){
    stop('Type Error : expt must be integer')
  }
  
  if(expt < 1){
    stop('Value Error : expt should be positive')
  }
  
  
  if(!(is.double(N)) ||  !(N %% 1 == 0)){
    stop('Type Error : N must be integer')
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
  
  stopifnot("input file not found" = file.exists(glue('testing/data/temp/{dataset}_expt_{expt}.rds')))
  
  message(glue('Running testing for Experiment {expt}'))
  # read simulated data for this experiment iteration and unload the list into variables
  zeallot::`%<-%`(c(sim_counts, sim_meta, effect_genes, null_genes, .), read_rds(glue('testing/data/temp/{dataset}_expt_{expt}.rds')))
  
  # get metadata for different scales of data representation
  zeallot::`%<-%`(c(cell_level_meta, bin_level_meta, sample_level_meta, cell_type_meta), prepare_metadata(sim_meta, N))
  clusters <- unique(cell_type_meta$cell_type)
  
  ## convert to pseudobulk counts
  binned_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_bin)
  celltype_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_cluster)
  ## assert that sample IDs match
  stopifnot(all(rownames(binned_counts) == bin_level_meta$sample_bin))
  stopifnot(all(rownames(celltype_counts) == cell_type_meta$sample_cluster))
  
  #===================================================================
  ## prepare normalized counts data for running argo
  #===================================================================
  # normalize counts using DESeq variance stabilizing transform
  # append centered, normalized counts to metadata object to use with fitting
  argo.input <- bin_level_meta %>%
    bind_cols(scale(t(varianceStabilizingTransformation(t(binned_counts))), center = T, scale = F))
  
  
  # list of different baseline testing methods for benchmarking
  diff.testing.methods <- list(
    DESeq2 = run_deseq, 
    edgeR = run_edger, 
    Limma = run_limma, 
    tradeSeq = run_tradeseq, 
    Lamian = run_lamian)
  
  #======================================================
  ## perform differential testing
  #======================================================
  ## capture the current environment to pass on to differential testing functions
  env <- environment() 
  # res.list <- list()
  ## call differential testing passing the local environment
  res.list <- map(diff.testing.methods, ~ .x(env))
  
  ## run differential testing with ARGO using RBF kernels
  res.list[['argo_shared_noise']] <- run_argo(env, bin_meta = argo.input, FIT_CONDITIONAL_NOISE = FALSE)
  res.list[['argo_conditional_noise']] <- run_argo(env, bin_meta = argo.input, FIT_CONDITIONAL_NOISE = TRUE)
  
  write_rds(res.list, glue('testing/output/{dataset}/Expt_{expt}_results.rds'))
  #===================================================================
  ## evaluate performance of different testing algorithms
  #===================================================================
  eval_metrics <- do.call(rbind, lapply(res.list, compute_eval_metrics,
                      effect_genes, null_genes, fdr_thresh = 0.05))

  return(eval_metrics)
  
}