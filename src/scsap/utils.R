prepare_metadata <- function(df){
  ## create pseudo-bulk covariates
  ## columns of interest in obs are : age, donor_id, tissue, cell_type, clust_annot, (binned) isodepth,
  ## cell type names cannot have space in between -- replace with '_'
  # df1 <-  df %>% 
  #   mutate(tissue = str_replace_all(tissue, c(' ' = '_', 
  #                                             '-' = '_', 
  #                                             '/' = '_'))) %>%
  #   mutate(
  #     bin_observed = factor(ntile(binned_isodepth, n = N)), 
  #     slide_bin = paste(batch, bin_observed, sep = '_'), # bin within individual slides
  #     condition_bin = paste(age, bin_observed, sep = '_'), 
  #     sample_cluster = factor(paste(batch, tissue, sep = '_')),
  #     condition_cluster = factor(paste(age, tissue, sep = '_')), 
  #     bin_tradeseq = as.numeric(factor(binned_isodepth)), 
  #     slide_tradeseq = factor(paste(batch, bin_tradeseq, sep = '_')),
  #     condition_tradeseq = factor(paste(age, bin_tradeseq, sep = '_'))
  #   ) %>% 
  #   mutate(
  #     condition_bin = factor(condition_bin, 
  #                            levels =  gtools::mixedsort(unique(condition_bin))),
  #     slide_bin = factor(slide_bin, 
  #                         levels = gtools::mixedsort(unique(slide_bin)))) %>% 
  #   group_by(bin_observed) %>%
  #   mutate(coord_x = mean(binned_isodepth)) %>%
  #   ungroup()
  
  df1 <-  obs %>%
    mutate(tissue = str_replace_all(tissue, c(' ' = '_',
                                              '-' = '_',
                                              '/' = '_'))) %>%
    mutate(
      bin = as.numeric(factor(binned_isodepth)),
      slide_bin = paste(batch, bin, sep = '_'), # bin within individual slides
      condition_bin = paste(age, bin, sep = '_'),
    ) %>% 
    dplyr::rename(coord_x = binned_isodepth) %>%
    mutate(
      condition_bin = factor(condition_bin,
                             levels =  gtools::mixedsort(unique(condition_bin))),
      slide_bin = factor(slide_bin,
                          levels = gtools::mixedsort(unique(slide_bin))))
  
  ## go from cell-level to pseudobulk bin-level metadata
  df2 <- df1 %>% dplyr::count(
    age, batch, donor_id,
    bin, coord_x, 
    slide_bin, condition_bin, name = 'n_cells'
    ) %>% filter(n_cells > 25) %>%
    arrange(slide_bin)
  
  ## sample x condition information for coregionalization matrix
  df3 <- df2 %>% select(age, donor_id, batch) %>% distinct()
  
  # df4 <- df1 %>% select(
  #   sample_cluster, condition_cluster,
  #   age, donor_id, batch, tissue
  #   ) %>% 
  #   distinct() %>% 
  #   arrange(sample_cluster)
  
  # same as in first step, except retain binned isodepth scale for each slide
  # df5 <- df1 %>% 
  #     dplyr::select(
  #     age, batch, donor_id,
  #     bin_tradeseq, slide_tradeseq, 
  #     condition_tradeseq, coord_x =  binned_isodepth 
  #   ) %>% distinct() %>% 
  #   arrange(slide_tradeseq)
  
  return(list(df1, df2, df3))
  # return(list(df1, df2, df3, df4, df5))
  
}

run_lamian <- function(env){
  message('running tests with Lamian')
  ## lamian requires normalized counts data as input. 
  # i'm using the official scater package to perform log normalization
  sce <- SingleCellExperiment(assays = list(counts = t(env$cell_level_counts)))
  sce <- scater::logNormCounts(sce)
  
  ## lamian requires an intercept term for XDE
  pt <- env$cell_level_meta$isodepth; names(pt) <- env$cell_level_meta$spot_id
  design <- model.matrix(~ age, data = env$sample_level_meta)
  rownames(design) <- env$sample_level_meta$batch
  testvar <- 3
  
  cell_anno <- data.frame(
    Cell = env$cell_level_meta$spot_id,
    Sample = as.character(env$cell_level_meta$batch)
  )
  
  # I'm using the chisq test to generate p-values instead of using permutations
  ## I also only care about a global differential testing result over the whole pseudotime axis
  res.lamian <- lamian_test(expr = as.matrix(logcounts(sce)),
                            cellanno = cell_anno,
                            pseudotime = pt,
                            design = as.matrix(design), 
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
  sce <- fitGAM(counts = t(as.matrix(env$tradeseq_counts)), 
                pseudotime = env$tradeseq_meta$coord_x, 
                cellWeights = matrix(1, nrow = nrow(env$tradeseq_meta), ncol=1),
                conditions = env$tradeseq_meta$age, 
                nknots = 6, 
                verbose = F, 
                parallel = F)
  
  ## run pairwise DGE test between conditions
  condRes <- conditionTest(sce, pairwise=T) %>% rownames_to_column('gene')
  ## extract results for the specific contrast of interest
  res.tradeseq <- condRes %>% arrange(desc(waldStat_conds2vs3)) %>% 
    mutate(Rank = row_number(), 
           neglogp = -log(pvalue_conds2vs3 + 1e-16),
           fdr = p.adjust(pvalue_conds2vs3, method = 'BH')) %>% as_tibble()
  return(res.tradeseq)
}

run_limma <- function(env){
  message('running tests with Limma')
  ## prepare DGE object and design matrix
  Y <- DGEList(counts = t(env$cluster_level_counts), 
               samples = env$cluster_meta)
  Y <- calcNormFactors(Y)
  
  ## design matrix for contrast testing 
  design <- model.matrix(~ 0 + condition_cluster, data = env$cluster_meta)
  colnames(design) <- gsub("^condition_cluster", "", colnames(design))
  
  ## apply voom-transformation to generate the precision weights
  v <- voom(Y, design = design, plot = F)
  # fit initial linear model
  fit <- lmFit(v, design)
  
  ## extract results with contrasts
  res.list <- list()
  for(ctx in env$clusters){
    contrast_formula <- glue("X90wk_{ctx} - X4wk_{ctx}")
    con <- makeContrasts(contrasts = contrast_formula, levels = design)
    ## fit contrast to the linear model
    fit_contrast <- contrasts.fit(fit, contrasts =con)
    ## run empirical bayes moderation
    fit_ebayes <- eBayes(fit_contrast)
    ## save results to list
    res.list[[ctx]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
  }
  
  res.limma <- summarize_dge_results(res.list, rownames(Y$counts), 
                                     pval = 'P.Value', lfc = 'logFC')[['dge']]
  return(res.limma)
}

run_edger <- function(env){
  message('running tests with edgeR')
  Y <- DGEList(counts = t(env$cluster_level_counts), 
               samples = env$cluster_meta)
  Y <- calcNormFactors(Y)
  
  ## design matrix for contrast testing 
  design <- model.matrix(~ 0 + condition_cluster, data = env$cluster_meta)
  colnames(design) <- gsub("^condition_cluster", "", colnames(design))
  
  ## run edgeR
  Y <- estimateDisp(Y, design)
  fit <- glmQLFit(Y, design)
  
  ## extract results with contrasts
  res.list <- list()
  for(ctx in env$clusters){
    contrast_formula <- glue("X90wk_{ctx} - X24wk_{ctx}")
    con <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    ## run the test and extract results
    qlf <- glmQLFTest(fit, contrast = con)
    ## save results to list
    res.list[[ctx]] <- topTags(qlf, n = Inf, sort.by = "none")$table
  }
  
  res.edger <- summarize_dge_results(res.list, genes = rownames(Y$counts), 
                                     pval = 'PValue', lfc = 'logFC')[['dge']]
  return(res.edger) 
}

run_deseq <- function(env){
  message('running tests with DESeq2')
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(env$cluster_level_counts),
                                        colData = env$cluster_meta,
                                        design = ~ 0 + condition_cluster)
  dds <- DESeq(dds)
  res <- results(dds)
  
  ## gather results across all bins
  res.list <- list()
  for(ctx in env$clusters){
    res.list[[ctx]] <- results(dds, contrast= c('condition_cluster', paste('X90wk', ctx, sep='_'),
                                                paste('X24wk', ctx, sep='_')))
  }
  res.deseq <- summarize_dge_results(res.list, genes = rownames(dds), 
                                     pval = 'pvalue', lfc = 'log2FoldChange')[['dge']]
  
  return(res.deseq)
}