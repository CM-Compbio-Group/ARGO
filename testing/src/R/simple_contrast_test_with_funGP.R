## run the statistical tests with funGP
## compute the distribution of the test statistic under null hypothesis 
## for the given contrast vector. return analytical p-value and (index of) max deviation

simple_contrast_test_with_funGP <- function(gene, par, env){
  #' this function is going to be called G times for every iteration of simulation run -- i will access the 
  #' envonment variables directly instead of redefining parent.env every time
  
  ## pick out variable and parameters
  y <- env$bin_meta %>% pull(gene)
  ## derive the design matrix given the contrast vector 
  cc <- generate_cc_matrix(env$CONTRAST_VECTOR)
  
  ## unload the GP parameters
  bw = par[1]
  sigma_f = par[2]
  sigma_n = par[3:length(par)]
  Sigma <- diag(sigma_n[env$noise_diag] ^ 2)
  
  ## compute covariance env$KERNELs
  if(env$KERNEL == 'RBF'){
    if(env$N == env$M){ # if testing on same points
      Kxx <- Ktt <- Ktx <- sigma_f^2 * exp(-0.5 * env$dsqx / (bw^2))
    } else{
      Kxx <- sigma_f^2 * exp(-0.5 * env$dsqx / (bw^2))
      Ktt <- sigma_f^2 * exp(-0.5 * env$dsqt / (bw^2))
      Ktx <- sigma_f^2 * exp(-0.5* env$dsqtx / (bw^2))
    }
    
  }else if(env$KERNEL == 'Matern5'){
    ## compute explicitly for each case. only two computations for each env$KERNEL
    Kxx <- sigma_f^2 * (1 + (sqrt(5)*env$r_x / bw) + 5*env$dsqx/(3*bw^2)) * exp(-(sqrt(5)*env$r_x / bw))
    Ktt <- sigma_f^2 * (1 + (sqrt(5)*env$r_t / bw) + 5*env$dsqt/(3*bw^2)) * exp(-(sqrt(5)*env$r_t / bw))
    Ktx <- sigma_f^2 * (1 + (sqrt(5)*env$r_tx / bw) + 5*env$dsqtx/(3*bw^2)) * exp(-(sqrt(5)*env$r_tx / bw))
  } else if(env$KERNEL == 'Matern3'){
    Kxx <- sigma_f^2 * (1 + sqrt(3)*env$r_x / bw) * exp(-sqrt(3)*env$r_x / bw)
    Ktt <- sigma_f^2 * (1 + sqrt(3)*env$r_t / bw) * exp(-sqrt(3)*env$r_t / bw)
    Ktx <- sigma_f^2 * (1 + sqrt(3)*env$r_tx / bw) * exp(-sqrt(3)*env$r_tx / bw)
  }else{
    stop(sprintf("Unsupported env$KERNEL type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", env$KERNEL))
  }
  
  ## derive the covariance blocks under the full distribution
  ## evaluate the covariance blocks
  kk <- length(env$CONTRAST_VECTOR)
  S11 <- kronecker(diag(kk), Ktt)
  S12 <- kronecker(t(env$W), Ktx)
  
  S22_H1 <- kronecker(env$B, Kxx) + kronecker(Sigma, diag(env$N)) + 1e-6*diag(length(y))
  U_f <- chol(S22_H1)
  S22_H1_inv <- chol2inv(U_f)
  
  ## derive the transformation term
  A = S12 %*% S22_H1_inv
  
  ## the test statistic
  mu_hat <- A %*% y
  fungp_test_statistic <- drop(kronecker(t(env$CONTRAST_VECTOR), diag(env$M)) %*% mu_hat)
  idx <- which.max(abs(fungp_test_statistic)) # index of max deviation in test statistic
  
  ## derive the data covariance under the null model
  W_null <- tcrossprod(cc, env$W)
  B_contrast <- crossprod(W_null)
  S22_H0 <- kronecker(B_contrast, Kxx) + kronecker(Sigma, diag(env$N)) 
  
  ## derive covariance of the test statistic
  sigma_U_H0 <- A %*% tcrossprod(S22_H0, A)
  ccx <- kronecker(t(env$CONTRAST_VECTOR), diag(env$M))
  sigma_statistic <- ccx %*% tcrossprod(sigma_U_H0, ccx)
  # Force symmetry by averaging with its transpose
  sigma_statistic <- (sigma_statistic + t(sigma_statistic)) / 2
  
  ## global analytical p-value against null hypothesis
  fgp_res <- list('statistic' = fungp_test_statistic,
                  'neglogp' =  compute_analytical_pvalue(sigma_statistic, fungp_test_statistic),
                  'deviation' = fungp_test_statistic[idx],
                  'index' = idx,
                  'coord_max' = env$x_test[idx])
  
  if(env$RETURN_CONF_BANDS){
    klt <- KLT_CI_bands(sigma_statistic, 
                        n_samples = env$n_monte_carlo_simulations, 
                        conf.level = env$conf_level_alpha)
    return(append(fgp_res, klt))  
  }
  
  else{
    return(fgp_res)
  }
  
}