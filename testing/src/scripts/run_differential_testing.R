run_differential_testing <- function(expt, dataset, N = 25, M = 100){
  #' load simulated counts for the iteration and run differential testing across all algorithms
  #' simulated counts are first converted into pseudo-bulk counts with desired number of bins along the isodepth
  #' raw / normalized pseudo-bulk counts are used as input for the different testing algorithms
  #' differential testing results are returned as a list to caller function
  #' 
  #' @param expt experiment iteration
  #' @param dataset dataset label 
  #' @param N number of observation points along pseudotime
  #' @param M number of test points along pseudotime for differential testing with Argo
  #' Expression counts etc are read from memory, not passed as arguments: 
  #'  a list of the generated counts, metadata, and the set of effect and null genes
  #' 
  #' @returns list of evaluation metrics across differential testing methods on this simulation run
  #' 
  
  message(glue('Running testing for Experiment {expt}'))
  # read simulated data for this experiment iteration and unload the list into variables
  c(sim_counts, sim_meta, effect_genes, null_genes) %<-% read_rds(glue('testing/data/temp/{dataset}_expt_{expt}.rds'))
  
  # get metadata for different scales of data representation
  c(cell_level_meta, bin_level_meta, sample_level_meta) %<-% prepare_metadata(sim_meta, N, M)
  
  ## convert to pseudobulk counts
  binned_counts <- rowsum(t(sim_counts), group = cell_level_meta$sample_bin)
  stopifnot(all(rownames(binned_counts) == bin_level_meta$sample_bin))
  
  # list of different baseline testing methods for benchmarking
  diff.testing.methods <- list(
    DESeq2 = run_deseq, 
    edgeR = run_edger, 
    Limma = run_limma, 
    tradeSeq = run_tradeseq, 
    Lamian = run_lamian)
  
  #======================================================
  ## define the objects necessary for differential testing
  #======================================================
  # load commonly used objects into this environment so I don't have to keep making them multiple times
  dds <- DESeqDataSetFromMatrix(countData = t(binned_counts),
                                colData = bin_level_meta,
                                design = ~ 0 + condition_bin)
  ## perform counts normalization
  vsd <- varianceStabilizingTransformation(dds, blind = F)
  genes <- rownames(vsd) # list of genes in this experiment
  
  # for edgeR and limma
  Y <- DGEList(counts = t(binned_counts), 
               samples = bin_level_meta)
  Y <- calcNormFactors(Y)
  
  ## design matrix for contrast testing 
  design <- model.matrix(~ 0 + condition_bin, data = bin_level_meta)
  colnames(design) <- gsub("^condition_bin", "", colnames(design))
  
  #======================================================
  ## perform differential testing
  #======================================================
  ## capture the current environment to pass on to differential testing functions
  env <- environment() 
  # res.list <- list()
  ## call differential testing passing the local environment
  res.list <- map(diff.testing.methods, ~ .x(env))
  
  ## run differential testing with ARGO using RBF kernels
  res.list[['argo_shared_noise']] <- run_argo(env, FIT_CONDITIONAL_NOISE = FALSE)
  res.list[['argo_conditional_noise']] <- run_argo(env, FIT_CONDITIONAL_NOISE = TRUE)
  
  write_rds(res.list, glue('testing/output/{dataset}/Expt_{expt}_results.rds'))
  #===================================================================
  ## evaluate performance of different testing algorithms
  #===================================================================
  eval_metrics <- do.call(rbind, lapply(res.list, compute_eval_metrics, 
                      effect_genes, null_genes, fdr_thresh = 0.05))
  
  return(eval_metrics)
  
}



compute_ranked_auprc <- function(ranked_list, signature){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  positive_score <- scores[binary_labels]
  negative_score <- scores[!binary_labels]
  pr_obj <- pr.curve(scores.class0 = positive_score,
                     scores.class1 = negative_score, 
                     curve = T)
  return(pr_obj)
}



compute_eval_metrics <- function(df, effect_genes, null_genes, fdr_thresh = 0.05){
  called <- filter(df, fdr < fdr_thresh) %>% pull(gene)
  not_called <- filter(df, fdr >= fdr_thresh) %>% pull(gene)
  TP = length(intersect(called, effect_genes))
  FP = length(intersect(called, null_genes))
  TN = length(intersect(not_called, null_genes))
  FN = length(intersect(not_called, effect_genes))
  
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  false_positive_rate = FP/ (FP + TN)
  false_negative_rate = FN/ (FN + TP)
  
  precision = TP / (TP + FP)
  empirical_FDR = FP / (TP + FP)
  
  accuracy = (TP + TN) /  (TP + FP + TN + FN)
  F1 = (2 * TP) /  (2 * TP + FP + FN)
  
  probj = compute_ranked_auprc(df %>% arrange(desc(neglogp)) %>% pull(gene), signature = effect_genes)
  AUPRC = round(probj$auc.integral, 4)
  
  return(c('TP' = TP, 'FP' = FP, 'TN' = TN, 'FN' = FN, 
                          'sensitivity' = sensitivity, 
                          'specificity' = specificity, 
                          'false_positive_rate' = false_positive_rate,
                          'false_negative_rate' = false_negative_rate,
                          'precision' = precision,
                          'empirical_FDR' = empirical_FDR,
                          'accuracy' = accuracy,
                          'F1' = F1, 
                          'AUPRC' = AUPRC)
  )
}



prepare_metadata <- function(df, N, M){
  ## create pseudo-bulk covariates
  df1 <-  df %>% 
    mutate(
      bin_observed = factor(ntile(pseudotime, n = N)), 
      bin_test = factor(ntile(pseudotime, n = M)), 
      sample_bin = paste(sampleID, bin_observed, sep = '_'),
      condition_bin = paste(condition, bin_observed, sep = '_')
    ) %>% 
    mutate(
      condition_bin = factor(condition_bin, 
                             levels =  gtools::mixedsort(unique(condition_bin))),
      sample_bin = factor(sample_bin, 
                          levels = gtools::mixedsort(unique(sample_bin)))) %>% 
    group_by(bin_observed) %>%
    mutate(coord_x = mean(pseudotime)) %>%
    group_by(bin_test) %>%
    mutate(test_x = mean(pseudotime)) %>%
    ungroup()
  
  
  ## go from cell-level to pseudobulk bin-level metadata
  df2 <- df1 %>% dplyr::select(
    condition, replicate, sampleID, 
    bin_observed, coord_x, 
    sample_bin, condition_bin 
  ) %>% 
    distinct() %>% 
    arrange(sample_bin)
  
  ## sample x condition information for coregionalization matrix
  df3 <- df2 %>% select(condition, sampleID) %>% distinct()
  
  return(list(df1, df2, df3))
}




run_lamian <- function(env){
  message('running tests with Lamian')
  ## lamian requires normalized counts data as input. 
  # i'll perform variance stabilizing transformation to normalize the data
  ## lamian requires an intercept term for XDE
  pt <- env$bin_level_meta$coord_x; names(pt) <- env$bin_level_meta$sample_bin
  design <- model.matrix(~ condition, data = env$sample_level_meta)
  rownames(design) <- env$sample_level_meta$sampleID
  testvar <- 2
  
  cell_anno <- data.frame(
    Cell = env$bin_level_meta$sample_bin,
    Sample = as.character(env$bin_level_meta$sampleID)
  )
  
  # I'm using the chisq test to generate p-values instead of using permutations
  ## I also only care about a global differential testing result over the whole pseudotime axis
  res.lamian <- lamian_test(expr = as.matrix(assay(env$vsd)),
                            cellanno = cell_anno,
                            pseudotime = pt,
                            design = design, 
                            test.type= 'variable', 
                            testvar = testvar, 
                            test.method = 'chisq',
                            permuiter = 1, 
                            overall.only = T,
                            ncores = 1) 
  
  
  ## get differential dynamic genes statistics
  res.lamian <- res.lamian$statistics %>% 
    rownames_to_column('gene') %>% 
    as_tibble() %>% 
    mutate(neglogp = -log(pval.chisq.overall)) %>% 
    dplyr::rename(fdr = fdr.chisq.overall) %>% 
    arrange(desc(neglogp)) %>% 
    mutate(Rank = row_number())
  
  return(res.lamian)
}

run_tradeseq <- function(env){
  message('runnig tests with tradeSeq')
  ## fit GAM and get an sce object
  sce <- fitGAM(counts = t(env$binned_counts), 
                pseudotime = env$bin_level_meta$coord_x, 
                cellWeights = rep(1, nrow(env$bin_level_meta)),
                conditions = env$bin_level_meta$condition, 
                nknots = 6, 
                verbose = F, 
                parallel = T)
  
  ## run DGE test between conditions
  condRes <- conditionTest(sce) %>% rownames_to_column('gene')
  res.tradeseq <- condRes %>% arrange(desc(waldStat)) %>% 
    mutate(Rank = row_number(), 
           neglogp = -log(pvalue + 1e-16),
           fdr = p.adjust(pvalue, method = 'BH')) %>% as_tibble()
  return(res.tradeseq)
}

run_limma <- function(env){
  message('running tests with Limma')
  ## apply voom-transformation to generate the precision weights
  v <- voom(env$Y, design = env$design, plot = T)
  # fit initial linear model
  fit <- lmFit(v, env$design)
  
  ## extract results with contrasts
  res.list <- vector("list", env$N)
  for(bin in seq_len(env$N)){
    contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
    con <- makeContrasts(contrasts = contrast_formula, levels = env$design)
    ## fit contrast to the linear model
    fit_contrast <- contrasts.fit(fit, contrasts =con)
    ## run empirical bayes moderation
    fit_ebayes <- eBayes(fit_contrast)
    ## save results to list
    res.list[[bin]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
  }
  
  res.limma <- summarize_dge_results(res.list, env$genes, 
                                     pval = 'P.Value', lfc = 'logFC')[['dge']]
  return(res.limma)
}

run_edger <- function(env){
  message('running tests with edgeR')
  ## run edgeR
  Y <- estimateDisp(env$Y, env$design)
  fit <- glmQLFit(Y, env$design)
  
  ## extract results with contrasts
  res.list <- list()
  for(bin in seq_len(env$N)){
    contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
    con <- makeContrasts(contrasts = contrast_formula, levels = env$design)
    
    ## run the test and extract results
    qlf <- glmQLFTest(fit, contrast = con)
    ## save results to list
    res.list[[bin]] <- topTags(qlf, n = Inf, sort.by = "none")$table
  }
  
  res.edger <- summarize_dge_results(res.list, genes = env$genes, 
                                     pval = 'PValue', lfc = 'logFC')[['dge']]
  return(res.edger) 
}

run_deseq <- function(env){
  message('running tests with DESeq2')
  dds <- DESeq(env$dds)
  res <- results(dds)
  
  ## gather results across all bins
  res.list <- vector("list", env$N)
  for(bin in seq_len(env$N)){
    res.list[[bin]] <- results(dds, contrast= c('condition_bin', paste('condition_B', bin, sep='_'),
                                                paste('condition_A', bin, sep='_')))
  }
  res.deseq <- summarize_dge_results(res.list, genes = env$genes, 
                                     pval = 'pvalue', lfc = 'log2FoldChange')[['dge']]
  return(res.deseq)
}


