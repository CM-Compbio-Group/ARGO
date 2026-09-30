fit_scDesign3 <- function(sce){
  #' Train scDesign3 to simulate the input dataset and returns list of fitted objects
  #' 
  #' Sequentially constructs a datset object, fits regression models for each gene, 
  #' copula model to learn gene-gene correlation patterns and extracts out the estimated parameters
  #' The learned parameters and correlations are used to simulate new ground truth data with 
  #' specific differential expression patterns
  
  
  message('fitting scDesign3 parameters to provided dataset')
  sce.data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1"
  )
  
  sce.marginal <- fit_marginal(
    data = sce.data,
    predictor = "gene",
    mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
    sigma_formula = "1", 
    family_use = "nb",
    n_cores = 20,
    usebam = FALSE
  )
  
  sce.copula <- fit_copula(
    sce = sce,
    assay_use = "counts",
    marginal_list = sce.marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 20,
    input_data = sce.data$dat,
  )
  
  
  sce.para <- extract_para(
    sce = sce,
    marginal_list = sce.marginal,
    n_cores = 20,
    family_use = "nb",
    new_covariate = NULL,
    data = sce.data$dat
  )
  
  
  return(list(sce = sce,
              data = sce.data,
              copula = sce.copula,
              para = sce.para))
  
}


