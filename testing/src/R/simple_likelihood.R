## run optimization with log-transformed parameter values
## transform back within the functions, return gradients wrt log-parameters
## P_obs is permutation matrix from full data to observed data

## define likelihood and gradient functions adaptively

computeNegLogLikelihoodGradient <- function(par, y, envir){
  
  ## derive the gradients wrt parameters
  gradients <- numeric(length=length(par))
  
  ## read current parameter estimates
  bw <- exp(par[1])
  sigma_f <- exp(par[2])
  sigma_n <- exp(par[3:length(par)])
  Sigma <- diag(sigma_n[envir$noise_diag]^2)
  
  ## compute the covariance envir$KERNEL for the icm process
  if(envir$KERNEL == 'RBF'){
    K_kernel <- sigma_f^2 * exp(-0.5 * envir$dsqx / (bw^2)) 
    dKbw <- K_kernel * envir$dsqx / (bw^3)
  } else if(envir$KERNEL == 'Matern5'){
    matern_exp <- sqrt(5)*envir$r_x / bw
    K_kernel <- sigma_f^2 * (1 + matern_exp + 5*envir$dsqx / (3*bw^2)) * exp(-matern_exp)
    dKbw <- sigma_f^2 / (3*bw) * exp(-matern_exp) * (5*envir$dsqx / (bw^2)) * (1 + matern_exp)
  }else if(envir$KERNEL == 'Matern3'){
    matern_exp <- sqrt(3)*envir$r_x / bw
    K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
    dKbw <- sigma_f^2 * exp(-matern_exp) * (3*envir$dsqx / (bw^3))
  }else {
    stop(sprintf("Unsupported envir$KERNEL type: '%s'. Choose 'RBF' or 'Matern5' or 'Matern3.", envir$KERNEL))
  }
  
  ## this is the full covariance envir$KERNEL
  K <- kronecker(envir$B, K_kernel) + kronecker(Sigma, diag(envir$N)) + 1e-6 * diag(length(y))
  
  ## compute terms for gradient
  U <- chol(K)
  K_inv <- chol2inv(U)
  alpha <- K_inv %*% y  
  
  P <- tcrossprod(alpha) - K_inv
  
  ## gradients wrt covariance envir$KERNEL parameters
  ## tr AB = sum(A * t(envir$B)) (sum of elementwise product)
  ## gradient wrt bandwidth
  gradients[1] <- -0.5 * sum(P * t(kronecker(envir$B, dKbw)))
  
  ## gradient for noise term
  dKsf <- (2 / sigma_f) * K_kernel
  gradients[2] <- -0.5 * sum(P*t(kronecker(envir$B, dKsf)))
  
  ## gradient wrt sigma_n terms
  diag_P <- diag(P)
  for(dd in seq_along(sigma_n)){
    snx <- sigma_n[dd]
    # Identify indices corresponding to group 'dd' in the full matrix
    slide_mask <- rep(envir$noise_diag == dd, each = envir$N)
    gradients[2 + dd] <- -snx * sum(diag_P[slide_mask]) + 
      2*envir$LAMBDA*snx # regulraization term
  }
  
  return(gradients * exp(par)) ## return log-gradients 
}

## the likelihood function
computeNegLogLikelihood <- function(par, y, envir){
  
  
  # extract current parameter values
  bw <- exp(par[1]) ## bandwidth of envir$KERNEL
  sigma_f <- exp(par[2]) ## scale of squared-exponential noise
  sigma_n <- exp(par[3:length(par)]) ## noise per condition
  Sigma <- diag(sigma_n[envir$noise_diag]^2) ## observation noise in data
  
  
  if(envir$KERNEL == 'RBF'){
    K_kernel <- sigma_f^2 * exp(-0.5 * envir$dsqx / (bw^2)) 
  } else if(envir$KERNEL == 'Matern5'){
    matern_exp <- sqrt(5)*envir$r_x / bw
    K_kernel <- sigma_f^2 * (1 + matern_exp + 5 * envir$dsqx / (3 * bw^2)) * exp(-matern_exp)
  }else if(envir$KERNEL == 'Matern3'){
    matern_exp <- sqrt(3)*envir$r_x / bw
    K_kernel = sigma_f^2 * (1 + matern_exp) * exp(-matern_exp)
  }else {
    stop(sprintf("Unsupported envir$KERNEL type: '%s'. Choose 'RBF', 'Matern5' or 'Matern3.", envir$KERNEL))
  }
  
  
  ## cov envir$KERNEL for full data
  ## add a numerical-safety jitter term to the covariance envir$KERNEL
  K <- kronecker(envir$B, K_kernel) + kronecker(Sigma, diag(envir$N)) +  1e-6 * diag(length(y))
  
  ## supporting constants
  U <- chol(K)
  K_inv <- chol2inv(U)
  log_det_K <- 2*sum(log(diag(U)))
  
  ## compute three parts of the log-likelihood function
  n_obs <- length(y)
  t1 <- 0.5* c(crossprod(y, K_inv) %*% y)
  t2 <- 0.5 * n_obs * log(2*pi)
  t3 <- 0.5*log_det_K
  reg <- envir$LAMBDA * sum(sigma_n^2) ## add regularization term
  
  return(t1 + t2 + t3 + reg) ## return the nll value
  
}


