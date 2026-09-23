#' run DGE testing using funGP on simulated counts data 
## use single variance term instead of learning condition-specific noise
## fit to fewer test points, don't need to use all 2100 pseudotime values

## fit parameters 

run_argo <- function(KERNEL = 'RBF', FIT_CONDITIONAL_NOISE = FALSE){
  ## define cov kernel exponents
  ## pairwise distances and dist-squared
  r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
  r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
  r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2
  
  ## update noise_diag to have the same term for all conditions
  ifelse(FIT_CONDITIONAL_NOISE, 
         noise_diag <- as.numeric(factor(sample_meta$condition)), # condition specific noise
         noise_diag <- rep(1, length(sample_meta$sample)) # all samples have same noise term
  ) 
  
  n_noise <- length(unique(noise_diag)) # no. of noise parameters
  
  ## define the coregionalization matrix
  W <- model.matrix(~0 + factor(sample_meta$condition))
  B <- tcrossprod(W)
  
  ## call the parameter estimation script
  ## run gp parameter inference
  ## pairwise distances and dist-squared
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
                               x = x_observed,
                               y = y,
                               N = N, 
                               r = r_x,
                               r_squared = dsqx,
                               noise_diag = noise_diag,
                               B = B,
                               lambda = 1, # regularization for noise nuggets
                               KERNEL = KERNEL, # fit specific covariance-kernel to data
                               method = 'L-BFGS-B', 
                               lower = c(min_log_bw, rep(-3, 1 + n_noise)),
                               upper = c(max_log_bw,  rep(3, 1 + n_noise)),
                               control = list(maxit = 1000, trace = 2))
    
    kernelParams[[gene]] <- c(bandwidth = exp(optimResult$par[1]), 
                              sigma_f = exp(optimResult$par[2]), 
                              sigma_n = exp(optimResult$par[3:(3 + n_noise - 1)]))
    
    
  }
  
  kernel_params <- do.call('rbind', kernelParams)
  colnames(kernel_params) <- c('bandwidth', 'sigma_f', paste('sigma_n', 1:n_noise, sep = ''))
  # write_rds(kernel_params, file = glue('output/simulations/argo_noise_{FIT_CONDITIONAL_NOISE}_params.rds'))
  
  #===================================================================
  ## run differential testing with funGP
  #===================================================================
  ## control variables
  KLT_CONF_BAND <- TRUE
  CONTRAST_VECTOR <- c(-1, 1)
  
  funGP_fits <- list()
  for(gene in genes){
    message(gene)
    funGP_fits[[gene]]  <- simple_contrast_test_with_funGP(CONTRAST_VECTOR = CONTRAST_VECTOR, 
                                                           gene = gene, par = kernel_params[gene, ])
  }
  
  #===================================================================
  ## gather results and visualize confidence bands
  #===================================================================
  #-----------------------------------------------------
  ## summarize experiment results
  #-----------------------------------------------------
  ## compute exceedence metrics
  exc.metrics <- do.call('rbind', lapply(funGP_fits, compute_exceedence_metrics, x_test))
  
  ## aggregate results from functional testing
  ## append gene conditions
  funGP.df <- funGP_fits %>%
    map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'gene') %>%
    bind_cols(exc.metrics) %>% 
    mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
    mutate(stats = neglogp * sign(deviation), 
           fdr = p.adjust(exp(-neglogp), method = 'BH'),
           net_auec = auec_up - auec_down,
           net_poe = poe_up - poe_down) %>%
    arrange(desc(neglogp)) %>%
    mutate(Rank = row_number())
  
  funGP.df

    
  print(funGP.df %>% filter(fdr < 0.05) %>% with(table(deviation > 0)))
  print(funGP.df %>% filter(fdr < 0.05 & max_deviation > 0) %>% with(table(gene %in% null_genes))) ## how many control genes are being called differential
  return(funGP.df)  
  
}

