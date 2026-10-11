## run optimization with log-transformed parameter values
## transform back within the functions, return gradients wrt log-parameters
## projection is permutation matrix from full data to observed data

## define likelihood and gradient functions adaptively

computeNegLogLikelihoodGradient <- function(par, argo_env){
  argo_env$.current_params <- par
  with(argo_env, {
    ## derive the gradients wrt parameters
    gradients <- numeric(length=length(.current_params))
    
    ## read current parameter estimates
    bw <- exp(.current_params[1])
    sigma_f <- exp(.current_params[2])
    sigma_d <- ifelse(FIT_DONOR, exp(.current_params[3]), 0) 
    sigma_n <- exp(.current_params[4:length(.current_params)])
    
    Sigma <- diag(sigma_n[noise_diag]^2)
    
    ## compute the covariance KERNEL for the icm process
    if(KERNEL == 'RBF'){
      K_kernel <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2)) 
      dKbw <- K_kernel * dsqx / (bw^3)
    } else if(KERNEL == 'Matern5'){
      matern_exp <- sqrt(5)*r_x / bw
      K_kernel <- sigma_f^2 * (1 + matern_exp + 5*dsqx / (3*bw^2)) * exp(-matern_exp)
      dKbw <- sigma_f^2 / (3*bw) * exp(-matern_exp) * (5*dsqx / (bw^2)) * (1 + matern_exp)
    }else if(KERNEL == 'Matern3'){
      matern_exp <- sqrt(3)*r_x / bw
      K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
      dKbw <- sigma_f^2 * exp(-matern_exp) * (3*dsqx / (bw^3))
    }
    
    ## this is the full covariance KERNEL
    K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) 
    # project to observed data
    K_obs <- projection %*% tcrossprod(K, projection) + 
             sigma_d^2 * donor_identity + 
             1e-6 * diag(length(feature_vec))
    
    ## compute terms for gradient
    U <- chol(K_obs)
    K_inv <- chol2inv(U)
    alpha <- K_inv %*% feature_vec  
    
    P <- tcrossprod(alpha) - K_inv
    P_full <- crossprod(projection, P) %*% projection
    
    ## gradients wrt covariance KERNEL parameters
    ## tr AB = sum(A * t(B)) (sum of elementwise product)
    ## gradient wrt bandwidth
    gradients[1] <- -0.5 * sum(P_full * t(kronecker(B, dKbw)))
    
    ## gradient for noise term
    dKsf <- (2 / sigma_f) * K_kernel
    gradients[2] <- -0.5 * sum(P_full*t(kronecker(B, dKsf)))
    
    ## gradient wrt sigma_d term
    gradients[3] <- ifelse(FIT_DONOR, -sigma_d * sum(P * donor_identity), 0)
    
    ## gradient wrt sigma_n terms
    diag_P_full <- diag(P_full)
    for(dd in seq_along(sigma_n)){
      snx <- sigma_n[dd]
      # Identify indices corresponding to group 'dd' in the full matrix
      slide_mask <- rep(noise_diag == dd, each = N)
      gradients[3 + dd] <- -snx * sum(diag_P_full[slide_mask]) + 
        2*LAMBDA*snx # regulraization term
    }
    return(gradients * exp(.current_params)) ## return log-gradients 
  })
}

## the likelihood function
computeNegLogLikelihood <- function(par, argo_env){
  argo_env$.current_params <- par
  with(argo_env, {
    
    # browser()
    # extract current parameter values
    bw <- exp(.current_params[1]) ## bandwidth of KERNEL
    sigma_f <- exp(.current_params[2]) ## scale of squared-exponential noise
    sigma_d <- ifelse(FIT_DONOR, exp(.current_params[3]), 0) ## donor-variation 
    sigma_n <- exp(.current_params[4:length(.current_params)]) ## noise terms for each timepoint
    
    # noise matrix
    Sigma <- diag(sigma_n[noise_diag]^2)
    
    if(KERNEL == 'RBF'){
      K_kernel <- sigma_f^2 * exp(-0.5 * dsqx / (bw^2)) 
    } else if(KERNEL == 'Matern5'){
      matern_exp <- sqrt(5)*r_x / bw
      K_kernel <- sigma_f^2 * (1 + matern_exp + 5 * dsqx / (3 * bw^2)) * exp(-matern_exp)
    }else if(KERNEL == 'Matern3'){
      matern_exp <- sqrt(3)*r_x / bw
      K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
    }
    
    
    ## cov KERNEL for full data
    ## add a numerical-safety jitter term to the covariance KERNEL
    K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) 
    
    ## project to observed data
    K <- projection %*% tcrossprod(K, projection) + 
      sigma_d^2 * donor_identity +  
      1e-6 * diag(length(feature_vec))
    
    ## supporting constants
    U <- chol(K)
    K_inv <- chol2inv(U)
    log_det_K <- 2*sum(log(diag(U)))
    
    ## compute three parts of the log-likelihood function
    N_obs <- length(feature_vec)
    t1 <- 0.5* c(crossprod(feature_vec, K_inv) %*% feature_vec)
    t2 <- 0.5 * N_obs * log(2*pi)
    t3 <- 0.5*log_det_K
    reg <- LAMBDA * sum(sigma_n^2) ## add regularization term
    
    return(t1 + t2 + t3 + reg) ## return the nll value
  })
}


