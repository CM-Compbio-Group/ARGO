## run optimization with log-transformed parameter values
## transform back within the functions, return gradients wrt log-parameters
## P_obs is projection matrix from full data to observed data

computeNegLogLikelihood <- function(par, x, y, P_obs, N, 
                                    r, noise_diag, slides_per_group,
                                    B, lambda ){
  # extract current parameter values
  bw <- exp(par[1]) ## bandwidth of kernel
  sigma <- exp(par[2]) ## scale of cov kernel noise
  sigma_n <- exp(par[3:length(par)]) ## noise terms for each timepoint
  
  ## compute the matern-5/2 covariance kernel
  matern_exp <- sqrt(5)*r / bw
  K_kernel <- sigma^2 * (1 + matern_exp + matern_exp^2 / 3) * exp(-matern_exp)
  
  Sigma <- diag(sigma_n[noise_diag]^2)
  
  ## cov kernel for observations
  ## add a numerical-safety jitter term to the covariance kernel
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) + 1e-6 * diag(N * sum(slides_per_group))
  
  ## transform to observed data
  K <- P_obs %*% tcrossprod(K, P_obs)
  
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

computeNegLogLikelihoodGradient <- function(par, x, y, P_obs, N, r,
                                            noise_diag, slides_per_group, B,
                                            lambda){
  ## read current parameter estimates
  bw <- exp(par[1])
  sigma <- exp(par[2])
  sigma_n <- exp(par[3:length(par)])
  
  ## derive the gradients wrt parameters
  gradients <- numeric(length=length(par))
  
  ## compute the covariance kernel for the icm process
  matern_exp <- sqrt(5)*r / bw
  K_kernel <- sigma^2 * (1 + matern_exp + matern_exp^2 / 3) * exp(-matern_exp)
  Sigma <- diag(sigma_n[noise_diag]^2)
  
  ## this is the full covariance kernel
  K <- kronecker(B, K_kernel) + kronecker(Sigma, diag(N)) + 1e-6 * diag(N * sum(slides_per_group))
  
  ## project to observed kernel
  K_obs <- P_obs %*% tcrossprod(K, P_obs)
  
  ## compute terms for gradient
  U <- chol(K_obs)
  K_inv <- chol2inv(U)
  alpha <- K_inv %*% y  
  
  P <- tcrossprod(alpha) - K_inv
  P_full <- crossprod(P_obs, P) %*% P_obs
  
  ## gradients wrt covariance kernel parameters
  ## tr AB = sum(A * t(B)) (sum of elementwise product)
  ## gradient wrt bandwidth
  # dKbw <- (K_kernel * matern_exp / bw) + (sigma^2 * exp(-matern_exp) * (-matern_exp / bw) * (1 + 2/3*matern_exp))
  ## simplify further
  dKbw <- sigma^2 / (3*bw) * exp(-matern_exp) * matern_exp^2 * (1 + matern_exp)
  gradients[1] <- -0.5 * sum(P_full * t(kronecker(B, dKbw)))
  
  
  ## gradient for noise term
  dKsf <- (2 / sigma) * K_kernel
  gradients[2] <- -0.5 * sum(P_full * t(kronecker(B, dKsf)))
  
  ## gradient wrt sigma_n terms
  diag_P_full <- diag(P_full)
  for(dd in seq_along(sigma_n)){
    sxn <- sigma_n[dd]
    
    # Identify indices corresponding to group 'dd' in the full matrix
    slide_mask <- rep(noise_diag == dd, each = N)
    gradients[2 + dd] <- -sxn * sum(diag_P_full[slide_mask]) + 
      2*lambda*sxn # regulraization term
  }
  
  return(gradients * exp(par)) ## return log-gradients 
}

