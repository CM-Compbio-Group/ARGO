## compute the distribution of the test statistic under null hypothesis 
## for the given contrast vector. return analytical p-value and (index of) max deviation

contrast_test_with_funGP <- function(feature, argo_env){
  argo_env$.feature <- feature
  with(argo_env, {
    ## derive the design matrix given the contrast vector 
    cc <- generate_cc_matrix(CONTRAST_VECTOR)
    params <- kernel_params[.feature, ]
    feature_vec <- counts[, .feature]
    
    ## unload the GP parameters
    bw = params[1]
    sigma_f = params[2]
    sigma_d = ifelse(FIT_DONOR, params[3], 0)
    sigma_n = params[4:length(params)]
    Sigma <- diag(sigma_n[noise_diag] ^ 2)
    
    ## compute covariance kernel
    if(KERNEL == 'RBF'){
      if(N == M){ # if testing on same points
        Kxx <- Ktt <- Ktx <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2))
      } else{
        Kxx <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2))
        Ktt <- sigma_f^2 * exp(-0.5 * dsqt / (bw^2))
        Ktx <- sigma_f^2 * exp(-0.5* dsqtx / (bw^2))
      }
    }else if(KERNEL == 'Matern3'){
      if(N == M){
        Kxx <- Ktt <- Ktx <- sigma_f^2 * (1 + sqrt(3)*r_x / bw) * exp(-sqrt(3)*r_x / bw)
      }else{
        Kxx <- sigma_f^2 * (1 + sqrt(3)*r_x / bw) * exp(-sqrt(3)*r_x / bw)
        Ktt <- sigma_f^2 * (1 + sqrt(3)*r_t / bw) * exp(-sqrt(3)*r_t / bw)
        Ktx <- sigma_f^2 * (1 + sqrt(3)*r_tx / bw) * exp(-sqrt(3)*r_tx / bw)
      }
    }else if(KERNEL == 'Matern5'){
      if(N == M){
        Kxx <- Ktt <- Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) +
                                            5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
      }else{
        ## compute explicitly for each case. only two computations for each KERNEL
        Kxx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
        Ktt <- sigma_f^2 * (1 + (sqrt(5)*r_t / bw) + 5*dsqt/(3*bw^2)) * exp(-(sqrt(5)*r_t / bw))
        Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_tx / bw) + 5*dsqtx/(3*bw^2)) * exp(-(sqrt(5)*r_tx / bw))
      }
    }
    
    ## derive the covariance blocks under the full distribution
    ## evaluate the covariance blocks
    kk <- length(CONTRAST_VECTOR)
    S11 <- kronecker(diag(kk), Ktt)
    S12_obs <- tcrossprod(kronecker(t(W), Ktx), projection)
    
    S22_H1 <- kronecker(B, Kxx) + kronecker(Sigma, diag(N)) 
    S22_H1_obs <- projection %*% tcrossprod(S22_H1, projection) + sigma_d^2*cov_donor
    U_f <- chol(S22_H1_obs)
    S22_H1_inv <- chol2inv(U_f)
    
    ## derive the transformation term
    A = S12_obs %*% S22_H1_inv
    
    ## the test statistic
    mu_hat <- A %*% feature_vec
    fungp_test_statistic <- drop(kronecker(t(CONTRAST_VECTOR), diag(M)) %*% mu_hat)
    index <- which.max(abs(fungp_test_statistic)) # index of max deviation in test statistic
    
    
    ## derive the data covariance under the null model
    W_null <- tcrossprod(cc, W)
    B_contrast <- crossprod(W_null)
    S22_H0 <- kronecker(B_contrast, Kxx) + kronecker(Sigma, diag(N)) 
    
    ## project into observed data
    S22_H0_obs <- projection %*% tcrossprod(S22_H0, projection) + sigma_d^2*cov_donor
    
    
    ## derive covariance of the test statistic
    sigma_U_H0 <- A %*% tcrossprod(S22_H0_obs, A)
    ccx <- kronecker(t(CONTRAST_VECTOR), diag(M))
    sigma_statistic <- ccx %*% tcrossprod(sigma_U_H0, ccx)
    # Force symmetry by averaging with its transpose
    sigma_statistic <- (sigma_statistic + t(sigma_statistic)) / 2
    
    ## global analytical p-value against null hypothesis
    fgp_res <- list('statistic' = fungp_test_statistic,
                    'neglogp' =  compute_analytical_pvalue(sigma_statistic, fungp_test_statistic),
                    'deviation' = fungp_test_statistic[index]
                    )
    
    if(RETURN_CONF_BANDS){
      klt <- KLT_CI_bands(sigma_statistic, 
                          n_samples = n_monte_carlo_simulations, 
                          conf.level = conf_level_alpha)
      return(append(fgp_res, klt))  
    }
    else{
      return(fgp_res)
    }
  })
}
   
