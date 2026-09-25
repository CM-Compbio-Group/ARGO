run_argo <- function(env, FIT_CONDITIONAL_NOISE = FALSE, KERNEL = 'RBF'){
  #' run differential testing with argo, first by fitting GP parameters for each gene and then 
  #' running differential testing. 
  #' 
  #' @param env environment of caller to directly access upstream variables
  #' @param FIT_CONDITIONAL_NOISE whether to fit likelihood with condition-specific noise terms or a shared nugget
  #' @param KERNEL what covariance kernel to use for fitting GPs to data
  #' @returns data frame with differential testing results for the genes
  
  message(glue("Running ARGO with conditional noise set to {FIT_CONDITIONAL_NOISE}"))
  
  # set parent environment
  rlang::env_poke_parent(environment(), new_env = env)
  
  
  ## read variables into the current environment (that i want to pass to nested functions)
  N <- env$N
  M <- env$M
  genes <- env$genes
  
  ## control variables
  LAMBDA <- 1 # regularization for noise terms in GP fits
  KLT_CONF_BAND <- TRUE
  CONTRAST_VECTOR <- c(-1, 1)
  RETURN_CONF_BANDS = TRUE
  conf_level_alpha = 0.95
  n_monte_carlo_simulations = 1000 # no.of monte carlo samples
  #===================================================================
  ## prepare normalized counts data for running argo
  #===================================================================
  # normalize counts using DESeq variance stabilizing transform
  # append centered, normalized counts to metadata object to use with fitting
  bin_meta <- bin_level_meta %>%
    bind_cols(scale(t(assay(vsd)), center = T, scale = F))
  
  #===================================================================
  ## prepare covariates for fitting and testing
  #===================================================================
  ## get variables for fitting continuous-time models
  x_observed <- sort(unique(cell_level_meta$coord_x)); 
  x_test <- sort(unique(cell_level_meta$test_x)); 
  
  ## define cov kernel exponents
  ## pairwise distances and dist-squared
  r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
  r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
  r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2
  
  #===================================================================
  ## prepare co-regionalization and noise matrices for GP fitting
  #===================================================================
  ## update noise_diag to have the same term for all conditions
  ifelse(FIT_CONDITIONAL_NOISE, 
         noise_diag <- as.numeric(factor(sample_level_meta$condition)), # condition specific noise
         noise_diag <- rep(1, length(sample_level_meta$sampleID)) # all samples have same noise term
  ) 
  n_noise <- length(unique(noise_diag)) # no. of noise parameters
  ## define the coregionalization matrix
  W <- model.matrix(~0 + factor(sample_level_meta$condition))
  B <- tcrossprod(W)
  #===========================================
  ## capture the variables in environment
  #===========================================
  argo_env <- environment()
  #===================================================================
  ## parameter estimation for each gene
  #===================================================================
  min_bw <- min(dist(x_observed)) * 0.5
  min_log_bw <- log(min_bw)
  max_bw <- 3*diff(range(x_observed))
  max_log_bw <- log(max_bw)
  
  kernelParams <- list()
  for(gene in genes){
    y <- bin_meta %>% pull(gene)
    log_bw <- runif(1, log(.1), log(.5))
    
    # compute empirical variation in expt1 by condition group
    emp_sd <- bin_meta %>% 
      {if(FIT_CONDITIONAL_NOISE) group_by(., condition) else .} %>%
      summarize(sd_val = sd(!!sym(gene))) %>% 
      pull(sd_val)
    
    log_sigma_f <-  log(mean(emp_sd) * sqrt(0.7)) ## 70% variance to signal
    log_sigma_n <- log(emp_sd * sqrt(0.3)) ## 30% variance to noise
    parInt <- c(log_bw, log_sigma_f, log_sigma_n)
    
    optimResult = stats::optim(par = parInt, 
                               fn = computeNegLogLikelihood, 
                               gr = computeNegLogLikelihoodGradient,
                               y = y, # pass the variable being fit
                               env = argo_env, # pass the current environment
                               method = 'L-BFGS-B', 
                               lower = c(min_log_bw, rep(-3, 1 + n_noise)),
                               upper = c(max_log_bw,  rep(3, 1 + n_noise)),
                               control = list(maxit = 1000, trace = 0)
                               )
    
    kernelParams[[gene]] <- c(bandwidth = exp(optimResult$par[1]), 
                              sigma_f = exp(optimResult$par[2]), 
                              sigma_n = exp(optimResult$par[3:(3 + n_noise - 1)]))
    
    
  }
  
  kernel_params <- do.call('rbind', kernelParams)
  colnames(kernel_params) <- c('bandwidth', 'sigma_f', paste('sigma_n', 1:n_noise, sep = ''))
  #===================================================================
  ## run differential testing with funGP
  #===================================================================
  funGP_fits <- list()
  for(gene in genes){
    funGP_fits[[gene]]  <- simple_contrast_test_with_funGP(gene = gene, 
                                                           par = kernel_params[gene, ], # pass fitted params for the gene
                                                           env = argo_env)
  }
  
  #===================================================================
  ## gather results and visualize confidence bands
  #===================================================================
  ## aggregate results from functional testing
  funGP.df <- funGP_fits %>%
    map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'gene') %>%
    mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
    mutate(stats = neglogp * sign(deviation), 
           fdr = p.adjust(exp(-neglogp), method = 'BH')
           ) %>%
    arrange(desc(neglogp)) %>%
    mutate(Rank = row_number())
  
  return(funGP.df)
  
}

