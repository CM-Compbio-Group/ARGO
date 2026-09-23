## run the statistical tests with funGP
## compute the distribution of the test statistic under null hypothesis 
## for the given contrast vector. return analytical p-value and (index of) max deviation

simple_contrast_test_with_funGP <- function(CONTRAST_VECTOR, gene, par,
                                           return_conf_band = TRUE,
                                           alpha = 0.95,
                                           n_mc = 1000 # no.of monte carlo samples
){
  
  
  
  ## pick out variable and parameters
  y <- bin_meta %>% pull(gene)
  
  ## derive the design matrix given the contrast vector 
  cc <- generate_cc_matrix(CONTRAST_VECTOR)
  
  ## unload the GP parameters
  bw = par[1]
  sigma_f = par[2]
  sigma_n = par[3:length(par)]
  Sigma <- diag(sigma_n[noise_diag] ^ 2)
  
  ## compute covariance KERNELs
  if(KERNEL == 'RBF'){
    Kxx <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2))
    Ktt <- sigma_f^2 * exp(-0.5 * dsqt / (bw^2))
    Ktx <- sigma_f^2 * exp(-0.5* dsqtx / (bw^2))
  }else if(KERNEL == 'Matern5'){
    ## compute explicitly for each case. only two computations for each KERNEL
    Kxx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
    Ktt <- sigma_f^2 * (1 + (sqrt(5)*r_t / bw) + 5*dsqt/(3*bw^2)) * exp(-(sqrt(5)*r_t / bw))
    Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_tx / bw) + 5*dsqtx/(3*bw^2)) * exp(-(sqrt(5)*r_tx / bw))
  } else if(KERNEL == 'Matern3'){
    Kxx <- sigma_f^2 * (1 + sqrt(3)*r_x / bw) * exp(-sqrt(3)*r_x / bw)
    Ktt <- sigma_f^2 * (1 + sqrt(3)*r_t / bw) * exp(-sqrt(3)*r_t / bw)
    Ktx <- sigma_f^2 * (1 + sqrt(3)*r_tx / bw) * exp(-sqrt(3)*r_tx / bw)
  }else{
    stop(sprintf("Unsupported KERNEL type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", KERNEL))
  }
  
  ## derive the covariance blocks under the full distribution
  ## evaluate the covariance blocks
  kk <- length(CONTRAST_VECTOR)
  S11 <- kronecker(diag(kk), Ktt)
  S12 <- kronecker(t(W), Ktx)
  
  S22_H1 <- kronecker(B, Kxx) + kronecker(Sigma, diag(N)) + 1e-6*diag(length(y))
  U_f <- chol(S22_H1)
  S22_H1_inv <- chol2inv(U_f)
  
  ## derive the transformation term
  A = S12 %*% S22_H1_inv
  
  ## the test statistic
  mu_hat <- A %*% y
  fungp_test_statistic <- drop(kronecker(t(CONTRAST_VECTOR), diag(M)) %*% mu_hat)
  idx <- which.max(abs(fungp_test_statistic)) # index of max deviation in test statistic
  
  
  ## derive the data covariance under the null model
  W_null <- tcrossprod(cc, W)
  B_contrast <- crossprod(W_null)
  S22_H0 <- kronecker(B_contrast, Kxx) + kronecker(Sigma, diag(N)) 
  
  ## derive covariance of the test statistic
  sigma_U_H0 <- A %*% tcrossprod(S22_H0, A)
  ccx <- kronecker(t(CONTRAST_VECTOR), diag(M))
  sigma_statistic <- ccx %*% tcrossprod(sigma_U_H0, ccx)
  # Force symmetry by averaging with its transpose
  sigma_statistic <- (sigma_statistic + t(sigma_statistic)) / 2
  
  ## global analytical p-value against null hypothesis
  fgp_res <- list('statistic' = fungp_test_statistic,
                  'neglogp' =  compute_analytical_pvalue(sigma_statistic, fungp_test_statistic),
                  'deviation' = fungp_test_statistic[idx],
                  'index' = idx)
  
  if(return_conf_band){
    klt <- KLT_CI_bands(sigma_statistic, n_samples = n_mc, conf.level = alpha)
    return(append(fgp_res, klt))  
  }
  
  else{
    return(fgp_res)
  }
  
}