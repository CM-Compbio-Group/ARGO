#' compute evaluation metrics based on ground truth data
compute_ranked_auprc <- function(ranked_list, signature){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  positive_score <- scores[binary_labels]
  negative_score <- scores[!binary_labels]
  pr_obj <- pr.curve(scores.class0 = positive_score, scores.class1 = negative_score, curve = T)
  return(pr_obj)
}


## pass list of true positives and true negatives -- compute eval metrics by effect sizes
compute_eval_metrics <- function(df, effect_genes, null_genes){
  
  ## subset to only this combination -- remove genes with other effect sizes
  df <- filter(df, gene %in% c(effect_genes, null_genes))
  called <- filter(df, fdr < 0.05) %>% pull(gene)
  not_called <- filter(df, fdr >= 0.05) %>% pull(gene)
  TP = length(intersect(called, effect_genes))
  FP = length(setdiff(called, effect_genes))
  TN = length(setdiff(not_called, effect_genes))
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
  AUPRC = round(probj$auc.integral, 3)
  
  return(list(metrics = c('TP' = TP, 'FP' = FP, 'TN' = TN, 'FN' = FN, 
           'sensitivity' = sensitivity, 
           'specificity' = specificity, 
           'false_positive_rate' = false_positive_rate,
           'false_negative_rate' = false_negative_rate,
           'precision' = precision,
           'empirical_FDR' = empirical_FDR,
           'accuracy' = accuracy,
           'F1' = F1, 
           'AUPRC' = AUPRC), 
           probj = probj)
  )
}
