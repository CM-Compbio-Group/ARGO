run_limma <- function(counts, meta){
  message('running tests with Limma')
  ## prepare DGE object and design matrix
  Y <- DGEList(counts = t(counts), 
               samples = meta)
  Y <- calcNormFactors(Y)
  
  ## design matrix for contrast testing 
  design <- model.matrix(~ 0 + condition_key, data = meta)
  colnames(design) <- gsub("^condition_key", "", colnames(design))
  
  ## apply voom-transformation to generate the precision weights
  v <- voom(Y, design = design, plot = F)
  # fit initial linear model
  fit <- lmFit(v, design)
  
  ## extract results with contrasts
  res.list <- list()
  for(ctx in levels(meta[['tissue']])){
    contrast_formula <- glue("X90wk___{ctx} - X24wk___{ctx}")
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

run_edger <- function(counts, meta){
  message('running tests with edgeR')
  Y <- DGEList(counts = t(counts), 
               samples = meta)
  Y <- calcNormFactors(Y)
  
  ## design matrix for contrast testing 
  design <- model.matrix(~ 0 + condition_key, data = meta)
  colnames(design) <- gsub("^condition_key", "", colnames(design))
  
  ## run edgeR
  Y <- estimateDisp(Y, design)
  fit <- glmQLFit(Y, design)
  
  ## extract results with contrasts
  res.list <- list()
  for(ctx in levels(meta[['tissue']])){
    contrast_formula <- glue("X90wk___{ctx} - X24wk___{ctx}")
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

run_deseq <- function(counts, meta){
  message('running tests with DESeq2')
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(counts),
                                        colData = meta,
                                        design = ~ 0 + condition_key)
  dds <- DESeq(dds)
  res <- results(dds)
  
  
  ## gather results across all bins
  res.list <- list()
  for(ctx in levels(meta[["tissue"]])){
    res.list[[ctx]] <- results(dds, contrast= c('condition_key', paste('X90wk', ctx, sep='___'),
                                                paste('X24wk', ctx, sep='___')))
  }
  res.deseq <- summarize_dge_results(res.list, genes = rownames(dds), 
                                     pval = 'pvalue', lfc = 'log2FoldChange')[['dge']]
  return(res.deseq)
}


run_lamian <- function(counts, meta){
  message('running tests with Lamian')
  
  sample_meta <- metadata %>% distinct(sample_id, donor_id, condition) %>% 
    mutate(condition = fct_relevel(condition, '24wk')) # make 24wk base condition
  
  ## lamian requires normalized counts data as input. 
  # i'm using the official scater package to perform log normalization
  sce <- SingleCellExperiment(assays = list(counts = t(counts)))
  sce <- scater::logNormCounts(sce)
  
  ## lamian requires an intercept term for XDE
  pt <- meta[['pseudotime']]; names(pt) <- meta[['cell_id']]
  design <- model.matrix(~ condition, data = sample_meta)
  rownames(design) <- sample_meta[['sample_id']]
  testvar <- 3
  
  cell_anno <- data.frame(
    Cell = meta[["cell_id"]],
    Sample = as.character(meta[['sample_id']])
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


run_tradeseq <- function(counts, meta){
  message('runnig tests with tradeSeq')
  seq_meta <- meta %>% mutate(bin = as.numeric(as.factor(pseudotime))) %>%
    mutate(grid_key = paste(sample_id, bin, sep = '___'), 
          condition_key = paste(condition, bin, sep = '___')) %>%
    mutate(grid_key = factor(grid_key, 
                             levels = gtools::mixedsort(unique(grid_key))),
           condition_key = factor(condition_key, 
                                  levels = gtools::mixedsort(unique(condition_key)))) %>%
    arrange(grid_key)
  
  seq_counts <- rowsum(counts, group = seq_meta$grid_key)
  seq_meta <- seq_meta %>%  select(-c(cell_id, tissue)) %>% 
    distinct() %>% arrange(grid_key)
  stopifnot('seq_counts ordering does not match with metadata' = all(rownames(seq_counts) == seq_meta$grid_key))
  
  ## fit GAM and get an sce object
  sce <- fitGAM(counts = t(as.matrix(seq_counts)), 
                pseudotime = seq_meta[["pseudotime"]], 
                cellWeights = matrix(1, nrow = nrow(seq_meta), ncol=1),
                conditions = seq_meta[['condition']], 
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

