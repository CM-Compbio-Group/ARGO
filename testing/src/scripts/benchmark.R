#' use scDesign3 to simulate ground truth single cell expression data along trajectories
#' that can be used to benchmark argo against baseline methods for differential testing.
#' We start by using scDesign to fit parameters for each gene in the dataset
#' We then perturb the mean parameters along the trajectory to create different kinds of effects, 
#' and check the capacity of algorithms to recover these deviances


dataset <- 'PANCREAS'
sce <- read_rds(glue('testing/data/{dataset}_sce.rds'))
print(sce) ## 1000 genes, 2087 cells

sce_fits <- fit_scDesign3(sce)

## iterate and generate perturbed datasets and store list of truths
N_EXPERIMENTS <- 10
EXPT_params <- lapply(N_EXPERIMENTS, generate_perturbed_data)

## call data processing for each iteration in parallel. Sequentially run different testing 
## algorithms and store results. 

