#' set up and run experiment to simulate gene expression data for two conditions using scDesign3
#' use the pancreas endocrinogenesis dataset to fit parameters for gene expression model along trajectory
#' use the fitted models to simulate ground truth data for two conditions, perturbing expression mean along the 
#' trajectory with set amplitude and bandwidth parameters
#' 
#' Store simulated expression data across experiments in output/ folder
#' 
#' DGP. 
#' we are given number of gnees to include, what percent have to be perturbed.
#' perturbation strength and bandwidth.  


message('Simulating perturbed gene expression counts data along single cell pseudotime trajectory')

suppressMessages({
source('testing/src/scripts/00_load_packages_and_functions.R')
})

dataset <- 'PANCREAS'
FIT_GENE_PARAMS <- FALSE

## The scenarios that are going to be generated
## make a map of parameter values and experiments
N_EXPERIMENTS <- 5 # no. of replicates for each parameter combination
perturb_amplitude <- c(0.1, 0.25, 0.5, 0.75, 1) # magnitude of LFC in perturbation
perturb_bandwidth <- c(0.01, 0.025, 0.05, 0.075, 0.1) # width of the perturbation
expt_df <- crossing(amplitude = perturb_amplitude,
                    bandwidth = perturb_bandwidth,
                    iter = factor(paste('iter', seq(N_EXPERIMENTS), sep = '-'))) %>% 
  mutate(Experiment = row_number())

print(glue('Size of the simulation grid is {nrow(expt_df)}'))
write_csv(expt_df, file = 'testing/data/temp/experiment_params.csv')

if(FIT_GENE_PARAMS){
  message(glue("Fitting scDesign3 parameters for {dataset} dataset"))
  sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))
  print(sce) ## 1000 genes, 2087 cells
  sce_fits <- fit_scDesign3(sce) # learn mean and variance parameters for all genes
  write_rds(sce_fits, glue('testing/output/{dataset}_sce_fits.rds'))
}




message("Using fitted parameters to generate perturbation data across TWO conditions")
## define control variables for perturbation experiment
MASTER_SEED <- 148
N_SAMPLES = 3 # number of randomly sampled replicates for each condition 
N_GENES = 500 # number of genes to use
perturb_frac = 0.2 # fraction of genes that must have perturbations
sce_fits <- read_rds(glue('testing/output/{dataset}_sce_fits.rds'))

set.seed(MASTER_SEED)
## generate the datasets in parallel and save to folder 'testing/data/temp/'
message('generating perturbation datasets')
base::lapply(seq_len(nrow(expt_df)), generate_perturbed_data, expt_df, sce_fits,
             N_GENES, perturb_frac, N_SAMPLES, dataset) |> futurize(seed = TRUE)


