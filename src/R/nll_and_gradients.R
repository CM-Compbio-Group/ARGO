## run optimization with log-transformed parameter values
## transform back within the functions, return gradients wrt log-parameters
## P_obs is permutation matrix from full data to observed data

## define likelihood and gradient functions adaptively

computeNegLogLikelihoodGradient <- function(par, x, y, P_obs, N, r, r_squared,
                                            noise_diag, slides_per_group,
                                            B, cov_donor, lambda, 
                                            kernel = 'RBF', fit_donor = F){
  
  ## derive the gradients wrt parameters
  gradients <- numeric(length=length(par))
  
  ## read current parameter estimates
  bw <- exp(par[1])
  sigma_f <- exp(par[2])
  sigma_d <- ifelse(fit_donor, exp(par[3]), 0) 
  sigma_n <- exp(par[4:length(par)])
  Sigma <- diag(sigma_n[noise_diag]^2)
  
  ## compute the covariance kernel for the icm process
  if(kernel == 'RBF'){
    K_kernel <- sigma_f^2 * exp(-0.5 * r_squared / (bw^2)) 
    dKbw <- K_kernel * r_squared / (bw^3)
  } else if(kernel == 'Matern5'){
    matern_exp <- sqrt(5)*r / bw
    K_kernel <- sigma_f^2 * (1 + matern_exp + 5*r_squared / (3*bw^2)) * exp(-matern_exp)
    dKbw <- sigma_f^2 / (3*bw) * exp(-matern_exp) * (5*r_squared / (bw^2)) * (1 + matern_exp)
  }else if(kernel == 'Matern3'){
    matern_exp <- sqrt(3)*r / bw
    K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
    dKbw <- sigma_f^2 * exp(-matern_exp) * (3*r_squared / (bw^3))
  }else {
    stop(sprintf("Unsupported kernel type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", kernel))
  }
  
  ## this is the full covariance kernel
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) 
  # project to observed data
  K_obs <- P_obs %*% tcrossprod(K, P_obs) + sigma_d^2 * cov_donor + 1e-6 * diag(length(y))
  
  ## compute terms for gradient
  U <- chol(K_obs)
  K_inv <- chol2inv(U)
  alpha <- K_inv %*% y  
  
  P <- tcrossprod(alpha) - K_inv
  P_full <- crossprod(P_obs, P) %*% P_obs
  
  ## gradients wrt covariance kernel parameters
  ## tr AB = sum(A * t(B)) (sum of elementwise product)
  ## gradient wrt bandwidth
  gradients[1] <- -0.5 * sum(P_full * t(kronecker(B, dKbw)))
  
  ## gradient for noise term
  dKsf <- (2 / sigma_f) * K_kernel
  gradients[2] <- -0.5 * sum(P_full*t(kronecker(B, dKsf)))
  
  ## gradient wrt sigma_d term
  gradients[3] <- ifelse(fit_donor, -sigma_d * sum(P * cov_donor), 0)
  
  ## gradient wrt sigma_n terms
  diag_P_full <- diag(P_full)
  for(dd in seq_along(sigma_n)){
    snx <- sigma_n[dd]
    # Identify indices corresponding to group 'dd' in the full matrix
    slide_mask <- rep(noise_diag == dd, each = N)
    gradients[3 + dd] <- -snx * sum(diag_P_full[slide_mask]) + 
      2*lambda*snx # regulraization term
  }
  
  return(gradients * exp(par)) ## return log-gradients 
}

## the likelihood function
computeNegLogLikelihood <- function(par, x, y, P_obs, N, r, r_squared,
                                    noise_diag, slides_per_group,
                                    B, cov_donor, lambda, 
                                    kernel = 'RBF', fit_donor = F ){
  
  
  # extract current parameter values
  bw <- exp(par[1]) ## bandwidth of kernel
  sigma_f <- exp(par[2]) ## scale of squared-exponential noise
  sigma_d <- ifelse(fit_donor, exp(par[3]), 0) ## donor-variation 
  sigma_n <- exp(par[4:length(par)]) ## noise terms for each timepoint
  Sigma <- diag(sigma_n[noise_diag]^2)
  
  
  if(kernel == 'RBF'){
    K_kernel <- sigma_f^2 * exp(-0.5 * r_squared / (bw^2)) 
  } else if(kernel == 'Matern5'){
    matern_exp <- sqrt(5)*r / bw
    K_kernel <- sigma_f^2 * (1 + matern_exp + 5 * r_squared / (3 * bw^2)) * exp(-matern_exp)
  }else if(kernel == 'Matern3'){
    matern_exp <- sqrt(3)*r / bw
    K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
  }else {
    stop(sprintf("Unsupported kernel type: '%s'. Choose 'RBF', 'Matern5' or 'Matern3.", kernel))
  }
  
  
  ## cov kernel for full data
  ## add a numerical-safety jitter term to the covariance kernel
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) 
  
  ## project to observed data
  K <- P_obs %*% tcrossprod(K, P_obs) + sigma_d^2 * cov_donor +  1e-6 * diag(length(y))
  
  ## supporting constants
  U <- chol(K)
  K_inv <- chol2inv(U)
  log_det_K <- 2*sum(log(diag(U)))
  
  ## compute three parts of the log-likelihood function
  N_obs <- length(y)
  t1 <- 0.5* c(crossprod(y, K_inv) %*% y)
  t2 <- 0.5 * N_obs * log(2*pi)
  t3 <- 0.5*log_det_K
  reg <- lambda * sum(sigma_n^2) ## add regularization term
  
  return(t1 + t2 + t3 + reg) ## return the nll value
  
}


