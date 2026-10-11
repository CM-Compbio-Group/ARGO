## function to gather all differential testing results and compute the evaluation metrics

eval.metrics <- vector("list", length = nrow(expt_df))
for(expt in seq_len(nrow(expt_df))){
  message(expt)
  
  zeallot::`%<-%`(c(., ., effect_genes, null_genes, .), read_rds(glue('testing/data/temp/{dataset}_expt_{expt}.rds')))
  res.list <- read_rds(glue('testing/output/{dataset}/Expt_{expt}_results.rds'))
  eval.metrics[[expt]] <- do.call(rbind, lapply(res.list, compute_eval_metrics, 
                                        effect_genes, null_genes, fdr_thresh = 0.05))
}

write_rds(eval.metrics, 'testing/output/diff_testing_eval_metrics.rds')

