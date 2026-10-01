#' set up and run experiment to simulate gene expression data for two conditions using scDesign3
#' use the pancreas endocrinogenesis dataset to fit parameters for gene expression model along trajectory
#' use the fitted models to simulate ground truth data for two conditions, perturbing expression mean along the 
#' trajectory with set amplitude and bandwidth parameters
#' 
#' Store simulated expression data across experiments in output/ folder


message('Simulating perturbed gene expression counts data along single cell pseudotime trajectory')

source('testing/src/scripts/00_load_packages_and_functions.R')
dataset <- 'PANCREAS'

## make a map of parameter values and experiments
N_EXPERIMENTS <- 3 # no. of replicates for each parameter combination
perturb_amplitude <- c(0.3, 0.6, 1) # magnitude of LFC in perturbation
perturb_bandwidth <- c(0.03, 0.06, 0.1) # width of the perturbation
expt_df <- crossing(amplitude = perturb_amplitude,
                    bandwidth = perturb_bandwidth,
                    iter = factor(paste('iter', seq(N_EXPERIMENTS), sep = '-'))) %>% 
  mutate(Experiment = row_number())
write_csv(expt_df, file = 'testing/data/temp/experiment_params.csv')

message(glue("Fitting scDesign3 parameters for {dataset} dataset"))
sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))
print(sce) ## 1000 genes, 2087 cells
sce_fits <- fit_scDesign3(sce) # learn mean and variance parameters for all genes
write_rds(sce_fits, glue('testing/output/{dataset}_sce_fits.rds'))



message("Using fitted parameters to generate perturbation data across TWO conditions")
## define control variables for perturbation experiment
N_SAMPLES = 3 # number of randomly sampled replicates for each condition 
N_GENES = 200 # number of genes to use
perturb_frac = 0.2 # fraction of genes that must have perturbations
sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))

## generate the datasets in parallel 
message('generating perturbation datasets')
base::lapply(seq_len(nrow(expt_df)), generate_perturbed_data, expt_df, sce_fits,
             N_GENES, perturb_frac, N_SAMPLES, dataset) |> futurize(seed = TRUE)


