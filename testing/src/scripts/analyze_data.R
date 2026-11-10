#' use simulated reads to perform differential testing between the two conditions using ARGO and other baseline methods
#' Compute performance metrics for different methods and write to file in 'testing/output/temp/'
#' ARGO has a parameter N that defines how many bins we split the data into along the isodepth axis

suppressMessages({
  source('testing/src/scripts/00_load_packages_and_functions.R')
})

message('setting up global variables:')
MASTER_SEED <- 148
dataset <- 'PANCREAS'
expt_df <- vroom::vroom('testing/data/temp/experiment_params.csv')
N <- 25

message('running differential testing')
plan(multisession, workers = 8)
set.seed(MASTER_SEED)
start_time <- Sys.time()
eval.metrics <- base::lapply(seq_len(nrow(expt_df)), 
                             run_differential_testing, N = N, dataset) |> futurize(seed = TRUE)
end_time <- Sys.time()
print(end_time - start_time)
write_rds(eval.metrics, glue('testing/output/diff_testing_eval_metrics.rds'))
plan(sequential)