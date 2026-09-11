## run optimization with log-transformed parameter values
## transform back within the functions, return gradients wrt log-parameters
## P_obs is permutation matrix from full data to observed data

## define likelihood and gradient functions adaptively

computeNegLogLikelihoodGradient <- function(par, x, y, N, r, r_squared,
                                            noise_diag,
                                            B, lambda, 
                                            kernel = 'RBF'){
  
  ## derive the gradients wrt parameters
  gradients <- numeric(length=length(par))
  
  ## read current parameter estimates
  bw <- exp(par[1])
  sigma_f <- exp(par[2])
  sigma_n <- exp(par[3:length(par)])
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
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) + 1e-6 * diag(length(y))
  
  ## compute terms for gradient
  U <- chol(K)
  K_inv <- chol2inv(U)
  alpha <- K_inv %*% y  
  
  P <- tcrossprod(alpha) - K_inv
  
  
  ## gradients wrt covariance kernel parameters
  ## tr AB = sum(A * t(B)) (sum of elementwise product)
  ## gradient wrt bandwidth
  gradients[1] <- -0.5 * sum(P * t(kronecker(B, dKbw)))
  
  ## gradient for noise term
  dKsf <- (2 / sigma_f) * K_kernel
  gradients[2] <- -0.5 * sum(P*t(kronecker(B, dKsf)))
  
  
  ## gradient wrt sigma_n terms
  diag_P <- diag(P)
  for(dd in seq_along(sigma_n)){
    snx <- sigma_n[dd]
    # Identify indices corresponding to group 'dd' in the full matrix
    slide_mask <- rep(noise_diag == dd, each = N)
    gradients[2 + dd] <- -snx * sum(diag_P[slide_mask]) + 
      2*lambda*snx # regulraization term
  }
  
  return(gradients * exp(par)) ## return log-gradients 
}

## the likelihood function
computeNegLogLikelihood <- function(par, x, y, N, r, r_squared,
                                    noise_diag,
                                    B, lambda, 
                                    kernel = 'RBF' ){
  
  
  # extract current parameter values
  bw <- exp(par[1]) ## bandwidth of kernel
  sigma_f <- exp(par[2]) ## scale of squared-exponential noise
  sigma_n <- exp(par[3:length(par)]) ## noise terms for each timepoint
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
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) +  1e-6 * diag(length(y))
  
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


