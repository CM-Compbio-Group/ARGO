#' use simulated reads to perform differential testing between the two conditions using ARGO and other baseline methods
#' Compute performance metrics for different methods and write to file

expt_df <- vroom('testing/data/temp/experiment_params.csv')
message('running differential testing')
plan(multisession, workers = 8)
start_time <- Sys.time()
eval.metrics <- base::lapply(seq_len(nrow(expt_df)), 
                             run_differential_testing, N = 50, dataset) |> futurize(seed = TRUE)
end_time <- Sys.time()
print(end_time - start_time)
write_rds(eval.metrics, glue('testing/output/diff_testing_eval_metrics.rds'))
plan(sequential)