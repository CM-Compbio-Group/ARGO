prepare_metadata <- function(df, N){
  ## create pseudo-bulk covariates
  ## cell type names cannot have space in between -- replace with '_'
  df1 <-  df %>% 
    mutate(cell_type = str_replace_all(cell_type, c(' ' = '_', 
                                                    '-' = '_'))) %>%
    mutate(
      bin_observed = factor(ntile(pseudotime, n = N)), 
      sample_bin = paste(sampleID, bin_observed, sep = '_'),
      condition_bin = paste(condition, bin_observed, sep = '_'), 
      sample_cluster = factor(paste(sampleID, cell_type, sep = '_')),
      condition_cluster = factor(paste(condition, cell_type, sep = '_'))
    ) %>% 
    mutate(
      condition_bin = factor(condition_bin, 
                             levels =  gtools::mixedsort(unique(condition_bin))),
      sample_bin = factor(sample_bin, 
                          levels = gtools::mixedsort(unique(sample_bin)))) %>% 
    group_by(bin_observed) %>%
    mutate(coord_x = mean(pseudotime)) %>%
    ungroup()
  
  ## go from cell-level to pseudobulk bin-level metadata
  df2 <- df1 %>% dplyr::select(
    condition, sampleID, 
    bin_observed, coord_x, 
    sample_bin, condition_bin 
  ) %>% 
    distinct() %>% 
    arrange(sample_bin)
  
  ## sample x condition information for coregionalization matrix
  df3 <- df2 %>% select(condition, sampleID) %>% distinct()
  
  df4 <- df1 %>% select(
    sample_cluster, condition_cluster,
    sampleID, condition, cell_type
    ) %>% 
    distinct() %>% 
    arrange(sample_cluster)
  
  return(list(df1, df2, df3, df4))
  
}


## compute overlap metrics between recovered intervals and the ground truth
compute_overlap_statistics <- function(intA, intB){
  if(is_empty(intA) && is_empty(intB)){
    stop('intervals to evaluate cannot be empty')
  }
  jaccard <- length(intersect(intA,intB)) / length(union(intA, intB))
  overlap.coeff <- length(intersect(intA,intB)) / min(length(intA), length(intB))
  dice.coeff <- 2 * length(intersect(intA,intB)) / (length(intA) + length(intB))
  precision <- length(intersect(intA,intB)) / length(intA)
  recall <- length(intersect(intA,intB)) / length(intB)
  return(list(jaccard_index = jaccard,
              overlap_coefficient = overlap.coeff,
              dice_coefficient = dice.coeff,
              precision = precision,
              recall = recall))
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
  if(!(is.data.frame(df))){
    stop('Type Error : df must be a data frame')
  }
  if(!(is.numeric(fdr_thresh))){
    stop('Type Error : fdr_thresh must be numeric')
  }
  if(fdr_thresh < 0 || fdr_thresh > 1){
    stop('Value Error : fdr_thresh must be between zero and one')
  }
  if(is_empty(effect_genes) || is_empty(null_genes)){
    stop('list of positive and negative genes cannot be zero')
  }
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
  
  # for computing PR curve, make sure genes are ordered by pvalues only, not by directional statistics
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



run_lamian <- function(env){
  message('running tests with Lamian')
  ## lamian requires normalized counts data as input. 
  # i'm using the official scater package to perform log normalization
  with(env, {
    sce <- SingleCellExperiment(assays = list(counts = sim_counts))
    sce <- scater::logNormCounts(sce)
    
    ## lamian requires an intercept term for XDE
    pt <- sim_meta$pseudotime; names(pt) <- sim_meta$cell_id
    design <- model.matrix(~ condition, data = sample_level_meta)
    rownames(design) <- sample_level_meta$sampleID
    testvar <- 2
    
    cell_anno <- data.frame(
      Cell = sim_meta$cell_id,
      Sample = as.character(sim_meta$sampleID)
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
    # print(res.lamian %>% filter(neglogp > -log(.05)) %>% with(table(gene %in% null_genes)))
    return(res.lamian)
  })
  
}

run_tradeseq <- function(env){
  message('runnig tests with tradeSeq')
  ## fit GAM and get an sce object
  with(env, {
    sce <- fitGAM(counts = sim_counts, 
                  pseudotime = sim_meta$pseudotime, 
                  cellWeights = matrix(1, nrow = nrow(sim_meta), ncol=1),
                  conditions = sim_meta$condition, 
                  nknots = 6, 
                  verbose = F, 
                  parallel = F)
    
    ## run DGE test between conditions
    condRes <- conditionTest(sce) %>% rownames_to_column('gene')
    res.tradeseq <- condRes %>% arrange(desc(waldStat)) %>% 
      mutate(Rank = row_number(), 
             neglogp = -log(pvalue + 1e-16),
             fdr = p.adjust(pvalue, method = 'BH')) %>% as_tibble()
    return(res.tradeseq)
  })
}
  
run_limma <- function(env){
    message('running tests with Limma')
  with(env, {
    ## prepare DGE object and design matrix
    Y <- DGEList(counts = t(celltype_counts), 
                 samples = cell_type_meta)
    Y <- calcNormFactors(Y)
    
    ## design matrix for contrast testing 
    design <- model.matrix(~ 0 + condition_cluster, data = cell_type_meta)
    colnames(design) <- gsub("^condition_cluster", "", colnames(design))
    
    ## apply voom-transformation to generate the precision weights
    v <- voom(Y, design = design, plot = F)
    # fit initial linear model
    fit <- lmFit(v, design)
    
    ## extract results with contrasts
    res.list <- list()
    for(ctx in clusters){
      contrast_formula <- glue("condition_B_{ctx} - condition_A_{ctx}")
      con <- makeContrasts(contrasts = contrast_formula, levels = design)
      ## fit contrast to the linear model
      fit_contrast <- contrasts.fit(fit, contrasts =con)
      ## run empirical bayes moderation
      fit_ebayes <- eBayes(fit_contrast)
      ## save results to list
      res.list[[ctx]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
    }
    
    p.vals <- do.call(cbind, lapply(res.list, function(res) res[genes, 'P.Value']))
    return(p.vals)
    # res.limma <- summarize_dge_results(res.list, rownames(Y$counts), 
    #                                    pval = 'P.Value', lfc = 'logFC')[['dge']]
    # return(res.limma)
  })
}

run_edger <- function(env){
  message('running tests with edgeR')
  with(env, {
    Y <- DGEList(counts = t(celltype_counts), 
                 samples = cell_type_meta)
    Y <- calcNormFactors(Y)
    
    ## design matrix for contrast testing 
    design <- model.matrix(~ 0 + condition_cluster, data = cell_type_meta)
    colnames(design) <- gsub("^condition_cluster", "", colnames(design))
    
    ## run edgeR
    Y <- estimateDisp(Y, design)
    fit <- glmQLFit(Y, design)
    
    ## extract results with contrasts
    res.list <- list()
    for(ctx in clusters){
      contrast_formula <- glue("condition_B_{ctx} - condition_A_{ctx}")
      con <- makeContrasts(contrasts = contrast_formula, levels = design)
      
      ## run the test and extract results
      qlf <- glmQLFTest(fit, contrast = con)
      ## save results to list
      res.list[[ctx]] <- topTags(qlf, n = Inf, sort.by = "none")$table
    }

    p.vals <- do.call(cbind, lapply(res.list, function(res) res[genes, 'PValue']))
    return(p.vals)
    # res.edger <- summarize_dge_results(res.list, genes = rownames(Y$counts), 
    #                                    pval = 'PValue', lfc = 'logFC')[['dge']]
    # return(res.edger) 
  })
}

run_deseq <- function(env){
  message('running tests with DESeq2')
  with(env, {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(celltype_counts),
                                          colData = cell_type_meta,
                                          design = ~ 0 + condition_cluster)
    dds <- DESeq(dds)
    res <- results(dds)
    
    ## gather results across all bins
    res.list <- list()
    for(ctx in clusters){
      res.list[[ctx]] <- results(dds, contrast= c('condition_cluster', paste('condition_B', ctx, sep='_'),
                                                  paste('condition_A', ctx, sep='_')))
    }
    
    p.vals <- do.call(cbind, lapply(res.list, function(res) res[genes, 'pvalue']))
    return(p.vals)
    
    # res.deseq <- summarize_dge_results(res.list, genes = rownames(dds), 
    #                                    pval = 'pvalue', lfc = 'log2FoldChange')[['dge']]
    # return(res.deseq)
  })
  
}


