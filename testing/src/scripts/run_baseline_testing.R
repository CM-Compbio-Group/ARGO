#' load the scDesign3-simulated counts and perform differential testing with various approaches
#' I need to convert the data into pseudobulk binned counts first, do library size normalization and then feed to the algorithms
#' 
## scripts to gather DGE testing results

# load('output/simulations/deseq_objects.rdata')

run_baseline_methods <- function(dds, binned_counts, bin_level_meta, sample_meta, N){

  #===================================================================
  ## run differential testing with deseq
  #===================================================================
  message('running tests with DESeq2')
  design(dds) <- ~ 0 + condition_bin
  dds <- DESeq(dds)
  res <- results(dds)
  print(summary(res, alpha = .05))
  
  ## gather results across all bins
  res.list <- vector("list", N)
  for(bin in seq_len(N)){
    res.list[[bin]] <- results(dds, contrast= c('condition_bin', paste('condition_B', bin, sep='_'),
                                                          paste('condition_A', bin, sep='_')))
  }
  
  
  head(res.list[[1]])
  res.deseq <- summarize_dge_results(res.list, genes = rownames(dds), 
                                     pval = 'pvalue', lfc = 'log2FoldChange')[['dge']]
  
  ## summarize findings
  p1 <- get_ranked_dge_plot(res.deseq, method = 'DESeq2', markers = effect_genes)
  # with(res.deseq, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
  print(res.deseq %>% filter(fdr < 0.05) %>% with(table(max.lfc > 0)))
  print(res.deseq %>% filter(fdr < 0.05) %>% with(table(gene %in% null_genes))) ## how many control genes are being called differential
  
  #===================================================================
  ## run differential testing with edgeR
  #===================================================================
  message('running tests with edgeR')
  y <- DGEList(counts = t(binned_counts), 
               samples = bin_level_meta)
  
  ## normalize the counts
  y <- calcNormFactors(y)
  ## design matrix for contrasts
  design <- model.matrix(~ 0 + condition_bin, data = y$samples)
  colnames(design) <- gsub("^condition_bin", "", colnames(design))
  
  ## run edgeR
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  ## extract results with contrasts
  res.list <- list()
  for(bin in seq_len(N)){
    contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
    con <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    ## run the test and extract results
    qlf <- glmQLFTest(fit, contrast = con)
    ## save results to list
    res.list[[bin]] <- topTags(qlf, n = Inf, sort.by = "none")$table
  }
  
  # names(res.list[[1]])
  res.edger <- summarize_dge_results(res.list, genes = rownames(dds), pval = 'PValue', lfc = 'logFC')[['dge']]
  
  p2 <- get_ranked_dge_plot(res.edger, method = 'edgeR', markers = effect_genes)
  # with(res.edger, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
  print(res.edger %>% filter(fdr < 0.05) %>% with(table(max.lfc > 0)))
  print(res.edger %>% filter(fdr < 0.05) %>% with(table(gene %in% null_genes))) ## how many control genes are being called differential
  
  #===================================================================
  ## run differential testing with limma
  #===================================================================
  message('running tests with Limma')
  ## apply voom-transformation to generate the precision weights
  v <- voom(y, design = design, plot = T)
  # fit initial linear model
  fit <- lmFit(v, design)
  
  ## extract results with contrasts
  res.list <- vector("list", N)
  for(bin in seq_len(N)){
    contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
    con <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    ## fit contrast to the linear model
    fit_contrast <- contrasts.fit(fit, contrasts =con)
    
    ## run empirical bayes moderation
    fit_ebayes <- eBayes(fit_contrast)
    
    ## save results to list
    res.list[[bin]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
  }
  
  
  res.limma <- summarize_dge_results(res.list, rownames(dds), pval = 'P.Value', lfc = 'logFC')[['dge']]
  
  (p3 <- get_ranked_dge_plot(res.limma, method = 'limma-voom', markers = effect_genes))
  # with(res.limma, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
  print(res.limma %>% filter(fdr < .05) %>% with(table(max.lfc > 0)))
  print(res.limma %>% filter(fdr < .05) %>% with(table(gene %in% null_genes)))
  
  #===================================================================
  ## run differential testing with tradeSeq
  #===================================================================
  message('runnig tests with tradeSeq')
  ## fit GAM and get an sce object
  sce <- fitGAM(counts = counts(dds, normalized = F), 
                pseudotime = bin_level_meta$coord_x, 
                cellWeights = rep(1, nrow(bin_level_meta)),
                conditions = bin_level_meta$condition, 
                nknots = 6, 
                verbose = T, 
                parallel = T)
  
  
  ## run DGE test between conditions
  condRes <- conditionTest(sce)
  mean(condRes$pvalue <= 0.05, na.rm = TRUE)
  
  ## order genes by significance
  condRes <- condRes %>% rownames_to_column('gene') %>% arrange(desc(waldStat))
  # plot(seq(G), condRes$waldStat)
  
  
  res.tradeseq <- condRes %>% arrange(desc(waldStat)) %>% 
    mutate(Rank = row_number(), 
           neglogp = -log(pvalue + 1e-16),
           fdr = p.adjust(pvalue, method = 'BH')) %>% as_tibble()
  
  (p4 <- get_ranked_dge_plot(res.tradeseq, method = 'tradeSeq', markers = effect_genes))
  print(res.tradeseq %>% filter(fdr < .05) %>% with(table(gene %in% null_genes)))

  #===================================================================
  ## run differential testing with Lamian
  #===================================================================
  message('running tests wtih Lamian')
  ## get normalized counts as input 
  all(bin_level_meta$sample_bin == colnames(vsd))
  pt <- bin_level_meta$coord_x; names(pt) <- bin_level_meta$sample_bin
  design <- model.matrix(~ condition, data = sample_meta)
  rownames(design) <- sample_meta$sample
  testvar <- 2
  
  cell_anno <- data.frame(
    Cell = bin_level_meta$sample_bin,
    Sample = as.character(bin_level_meta$sample)
  )
  
  # get results for contrast controlling for donor effects
  res.lamian <- lamian_test(expr = as.matrix(assay(vsd)),
                     cellanno = cell_anno,
                     pseudotime = pt,
                     design = design, 
                     test.type= 'variable', 
                     testvar = testvar, 
                     test.method = 'chisq',
                     permuiter = 100, 
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
  
  (p5 <- get_ranked_dge_plot(res.lamian, method = 'Lamian', markers = effect_genes))
  print(res.lamian %>% filter(fdr < .05))
  print(res.lamian %>% filter(fdr < .05) %>% with(table(gene %in% null_genes)))
  

  ## visualize diff testing results across all methods
  pdf('figs/simx_{iter}.pdf', width = 22, height=6)
  plot(wrap_plots(p1 | p2 | p3 | p4 | p5))
  dev.off()
  
  return(list('deseq' = res.deseq,
              'edgeR' = res.edger,
              'limma' = res.limma,
              'tradeseq' = res.tradeseq,
              'lamian' = res.lamian)
  )
}