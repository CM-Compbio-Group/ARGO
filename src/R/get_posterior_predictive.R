## return posterior predictive plots for the genes given the parameter fits
get_posterior_predictive <- function(feature, argo_env){
  argo_env$.feature <- feature
  
  with(argo_env, {
    # browser()
    feature_vec <- counts[, .feature]
    params <- kernel_params[.feature, ]
    bw = params[1]
    sigma_f = params[2]
    sigma_d = ifelse(FIT_DONOR, params[3], 0)
    sigma_n = params[4:length(params)]
    
    n_slides <- length(levels(metadata$sample_id))
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
        Kxx <- Ktt <- Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
      }else{
      ## compute explicitly for each case. only two computations for each KERNEL
      Kxx <- sigma_f^2 * (1 + (sqrt(5)*r_x / bw) + 5*dsqx/(3*bw^2)) * exp(-(sqrt(5)*r_x / bw))
      Ktt <- sigma_f^2 * (1 + (sqrt(5)*r_t / bw) + 5*dsqt/(3*bw^2)) * exp(-(sqrt(5)*r_t / bw))
      Ktx <- sigma_f^2 * (1 + (sqrt(5)*r_tx / bw) + 5*dsqtx/(3*bw^2)) * exp(-(sqrt(5)*r_tx / bw))
      }
    }
    
    Kf <- kronecker(B, Kxx) + kronecker(Sigma, diag(N)) 
    Kf <- projection %*% tcrossprod(Kf, projection) + sigma_d^2*cov_donor
    U_f <- chol(Kf)
    Kf_inv <- chol2inv(U_f)
    
    ##--------------------------------------------------
    ## posterior mean for observations (UW^T | Y)
    ##--------------------------------------------------
    ## if we want donor-specific shifted curves for each mouse
    if(FIT_DONOR){
      donor_levels <- levels(metadata$sample_id)
      test_donors <- rep(metadata %>% distinct(sample_id) %>% pull(), each=M)
      Z_test <- model.matrix(~0 + factor(test_donors, levels = donor_levels))
      Z_obs <- model.matrix(~0 + metadata$sample_id)
      cov_donor_star <- tcrossprod(Z_test, Z_obs)
      Ktx_obs <- tcrossprod(kronecker(B, Ktx) , projection) + sigma_d^2*cov_donor_star
    }else{
      Ktx_obs <- tcrossprod(kronecker(B, Ktx) , projection)
    }
    ##--------------------------------------------------
    ## posterior predictives
    y_pred <- Ktx_obs %*% Kf_inv %*% feature_vec
    
    return(y_pred) # return posterior predictive counts
  })
}
