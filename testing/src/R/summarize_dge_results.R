## take list of differential expression results and return dataframe

summarize_dge_results <- function(res.list, genes, pval = 'p_val', lfc = 'avg_log2FC'){
  
  ## extract p-values and LFC from each test
  p.vals <- do.call(cbind, lapply(res.list, function(df) df[genes, pval]))
  lfc <- do.call(cbind, lapply(res.list, function(df) df[genes, lfc]))
  rownames(p.vals) <- rownames(lfc) <- genes
  
  ## get cauchy-combined neglog-pvals for each gene across tests
  nlpcx <- apply(p.vals, 1, function(p) cct_p.vals(p[!is.na(p)]))
  fdr = p.adjust(exp(-nlpcx), method = 'BH') 
 
  ## assign direction based on most significant deviation across regions
  max.idx <- apply(p.vals, 1, which.min)
  max.delta <- lfc[cbind(1:nrow(lfc), max.idx)] # deviation corresponding to most significant p-value
  names(max.delta) <- rownames(lfc)
  
  ## make separate summary dataframe
  dge <- data.frame(neglogp = nlpcx, 
                    fdr = fdr,
                    max.lfc = max.delta,
                    index = max.idx) %>% 
    mutate(stats = neglogp * sign(max.lfc)) %>% 
    arrange(desc(stats)) %>% 
    rownames_to_column('gene') %>% 
    mutate(Rank = row_number()) %>%
    as_tibble()
  
  return(list('dge' = dge,
              'p.vals' = p.vals,
              'lfc' = lfc))
}
