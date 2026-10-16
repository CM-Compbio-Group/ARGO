#' script to run ARGO -- fit parameters, run differential testing, save results and perform visualization
#' I want everything to be self-contained here. 
#' @param counts expression counts data
#' @param metadata metadata for each cell in dataset
#' @param features_to_fit list of features / cell types to test
#' @param N number of data points for fitting GP parameters
#' @param isodepth_ co-ordinates of observed data
#' @param test_coords coords to impute fits to
#' @param W design matrix for experiment
#' @param B coregionalization matrix
#' @param projection projection matrix into observed data
#' @param noise_diag ids for fitting condition-specific noise per condition
#' @param donor_identity donor identity matrix
#' @param KERNEL what covariance kernel to fit 
#' @param FIT_DONOR whether to fit donor level random intercepts
#' @param FIT_CONDITIONAL_NOISE = TRUE # logical -- whether to fit condition specific noise for GPs
#' @param LAMBDA <- 1 # regularization for noise terms in GP fits
#' @param CONTRAST_VECTOR <- c(0, -1, 1)
#' @param CONTRAST contrast being tested in words
#' @param RETURN_CONF_BANDS = TRUE # compute conf bands only if I am interested in recovering intervals
#' @param conf_level_alpha = 0.95
#' @param n_monte_carlo_simulations = 1000 # no.of monte carlo samples
#' @param kernel_params GP parameters for features to fit
#' 
#' @returns list of fitted kernel parameters and ARGO differential testing outputs


fit_ARGO <- function(data_env){
  orig_parent <- rlang::env_parent(environment())
  rlang::env_poke_parent(data_env, orig_parent) # make parent to the passed data environment the global environment
  rlang::env_poke_parent(environment(), new_env = data_env) # make data environment parent to the function environment
  
  ## ensure inputs are valid
  if(!is.double(N) && !(N %% 1 == 0)){
    stop('Type Error : N must be integer')
  }
  if(N < 10 || N > 100){
    stop('Value Error : number of intervals not in valid range')
  }
  if(!is.numeric(counts)){
    stop('Type Error : counts matrix must be numeric')
  }
  if(!(all(c('sample_id', 'bin', 'condition') %in% colnames(metadata)))){
    stop('Value Error : metadata must have columns "sample_id", "bin" and "condition')
  }
  stopifnot("counts matrix and metadata must be ordered" = all(rownames(counts) == metadata$grid_key))
  
  if(!is.numeric(CONTRAST_VECTOR)){
    stop('Type Error: contrast vector must be numeric vector')
  }
  if(sum(CONTRAST_VECTOR) != 0){
    stop('Value Error : contrast vector must sum to zero' )
  }
  if(!(1 %in% CONTRAST_VECTOR)){
    stop('Value Error : contrast vector must have a focal group' )
  }
  if(!is.numeric(isodepth_)){
    stop('Type Error: isodepth_ must be numeric vector')
  }
  if(!(is.character(KERNEL))){
    stop('Type Error: KERNEL must be a string')
  }
  if(!(KERNEL %in% c('RBF', 'Matern3', 'Matern5'))){
    stop('Value Error : KERNEL must be one of RBF, Matern3 or Matern5')
  }
  if(!is_logical(RETURN_CONF_BANDS)){
    stop('Type Error: RETURN_CONF_BANDS must be boolean')
  }
  if(!is_logical(FIT_DONOR)){
    stop('Type Error: FIT_DONOR must be boolean')
  }
  if(!is_logical(FIT_CONDITIONAL_NOISE)){
    stop('Type Error: FIT_CONDITIONAL_NOISE must be boolean')
  }
  if(conf_level_alpha < 0 || conf_level_alpha > 1){
    stop('Value Error: conf_level_alpha must be between 0 and 1')
  }
  if(!(is.numeric(projection) && is.matrix(projection))){
    stop('Type Error: projection must be a numeric matrix')
  }
  if(!(length(noise_diag) == length(unique(metadata$sample_id)))){
    stop('Value Error: noise vector does not align with metadata')
  }
  if(!is.double(n_monte_carlo_simulations) && !(n_monte_carlo_simulations %% 1 == 0)){
    stop('Type Error : n_monte_carlo_simulations must be integer')
  }
  if(!(is.character(features_to_fit))){
    stop('Type Error : list of features to fit must be strings')
  }
  if(!(all(features_to_fit %in% colnames(counts)))){
    stop('Value Error: desired features not found in counts matrix')
  }
  if(!(is.null(test_coords) || is_double(test_coords))){
    stop('Type Error: test_coords must be either NULL or a numeric vector')
  }
  if(!(is.numeric(kernel_params))){
    stop('Type Error : kernel params should be a numeric vector or matrix')
  }
  if(!(is.null(kernel_params)) && !(all(features_to_fit %in% rownames(kernel_params)))){
    stop('Value Error: GP parameters not available for all requested features')
  }

  
  
  ## define a new function environment 
  argo_env <- new.env(parent = data_env) # all further variables are now defined in this new environment
  
  with(argo_env, {
    ## set variables for fitting covariance kernels
    r_x <- as.matrix(dist(data_env$isodepth_));  dsqx <- r_x^2 # pairwise distances for fitting covariance kernels
    if(!(is.null(data_env$test_coords))){
      M <- length(data_env$test_coords)
      r_t <- as.matrix(dist(data_env$test_coords)); dsqt <- r_t^2
      r_tx <- rdist::cdist(data_env$test_coords, data_env$isodepth_); dsqtx <- r_tx^2
    }else{
      M <- data_env$N
      test_coords <- data_env$isodepth_
      r_t <- r_x; dsqt <- dsqx
      r_tx <- r_x; dsqtx <- dsqx
    }
    
    if(!data_env$FIT_CONDITIONAL_NOISE){
      noise_diag <- rep(1, length(levels(data_env$metadata$sample_id)))
    }else{
      noise_diag <- data_env$noise_diag
    }
  })
  
  if(is.null(kernel_params)){
    # only fit likelihood model when parameters are not available
    #--------------------------------
    ## set constraints for kernel bandwidth
    #--------------------------------
    min_bw <- min(dist(isodepth_)) * 0.5
    min_log_bw <- log(min_bw)
    max_bw <- 3*diff(range(isodepth_))
    max_log_bw <- log(max_bw)
    n_noise <- length(unique(argo_env$noise_diag)) # number of noise parameters to fit
    
    #===================================================================
    ## parameter estimation for each feature
    #===================================================================
    message(glue('now fitting gp parameters for {length(data_env$features_to_fit)} features'))
    
    kernel_params <- vector('list', length = length(data_env$features_to_fit))
    for(feature in data_env$features_to_fit){
      argo_env$feature_vec <- data_env$counts[, feature] # extract counts data for this feature
      
      # randomly initialize likelihood parameters 
      log_bw <- runif(1, log(.1), log(.5))
      # compute empirical variation in expt1 by condition group
      if(data_env$FIT_CONDITIONAL_NOISE){
        emp_sd <- tapply(argo_env$feature_vec, data_env$metadata$condition, sd, na.rm = T)
      }else{
        emp_sd <- sd(argo_env$feature_vec, na.rm=T)
      }
      
      log_sigma_f <-  log(mean(emp_sd) * sqrt(0.7)) ## 70% variance to signal
      log_sigma_n <- log(emp_sd * sqrt(0.3)) ## 30% variance to noise
      
      if(data_env$FIT_DONOR){
        sd_val <- tapply(argo_env$feature_vec, data_env$metadata$condition_key, sd, na.rm=T); # SD across donors within each bin
        log_sigma_d <- mean(sd_val, na.rm=T) # mean SD across all age groups
      }else{
        log_sigma_d <- 0
      }
      
      parInt <- c(log_bw, log_sigma_f, log_sigma_d, log_sigma_n)
      
      optimResult = stats::optim(par = parInt, 
                                 fn = computeNegLogLikelihood, 
                                 gr = computeNegLogLikelihoodGradient,
                                 argo_env = argo_env, # pass the current environment
                                 method = 'L-BFGS-B', 
                                 lower = c(min_log_bw, rep(-3, 2 + n_noise)),
                                 upper = c(max_log_bw,  rep(3, 2 + n_noise)),
                                 control = list(maxit = 1000, 
                                                trace = 0)
      )
      
      kernel_params[[feature]] <- c(bandwidth = exp(optimResult$par[1]), 
                                    sigma_f = exp(optimResult$par[2]), 
                                    sigma_d = exp(optimResult$par[3]), 
                                    sigma_n = exp(optimResult$par[4:(4 + n_noise - 1)]))
      
      
    }
    
    argo_env$kernel_params <- do.call('rbind', kernel_params)
    colnames(argo_env$kernel_params) <- c('bandwidth', 'sigma_f','sigma_d', paste('sigma_n', 1:n_noise, sep = ''))
    
  }else{
    message('GP kernel parameters provided by user')
    argo_env$kernel_params <- kernel_params
  }
    
  # visualize the posterior fits for genes
  message('now plotting posterior fits')
  posteriors <- lapply(data_env$features_to_fit, posterior_predictive_plots, argo_env)
  ggsave(wrap_plots(posteriors, ncol = 3) + 
         plot_layout(guides = "collect"), 
         filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_posterior_fits.pdf')),
         width=12,height=6)
  
  
  message('now running contrast tests')
  funGP_fits <- list()
  for(feature in data_env$features_to_fit){
    funGP_fits[[feature]]  <- contrast_test_with_funGP(feature, 
                                     argo_env)
  }
  
  
  
  message('now plotting confidence bands')
  ci_plots <- list()
  for(feature in data_env$features_to_fit){
    ci_plots[[feature]] <- plot_CI_bands(feature, funGP_fits[[feature]], argo_env)
  }
  ggsave( wrap_plots(ci_plots, ncol = 3) + 
            plot_layout(guides = "collect"), 
          filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_conf_bands.pdf')),
          width=10,height=5)
  
  
  if(is.null(data_env$kernel_params)){ # only return results if parameters are not known already
    return(list(resObjects = funGP_fits,
                params = argo_env$kernel_params)
    )  
  }  
  
  
}

