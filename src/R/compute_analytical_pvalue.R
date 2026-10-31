compute_analytical_pvalue <- function(sigma_statistic, stat, var_explained = 0.99, ridge_penalty = 1e-4){
  
  eigen_decomp <- eigen(sigma_statistic, symmetric = TRUE)
  evals <- eigen_decomp$values  
  
  ## restrict to number of eigenvalues explaining specified fraction of total variance
  cum_var <- cumsum(evals) / sum(evals)
  K <- sum(cum_var <= var_explained)
  
  # Ensure at least 1 dimension is kept if the first eigenvalue dominates >99%
  if (K == 0) K <- 1
  
  # threshold <- 1e-6*evals[1]
  # K <- sum(evals >= threshold)
  lambda_vals <- evals[1:K]
  phi <- eigen_decomp$vectors[, 1:K, drop = FALSE]
  
  ## project test statistic vector into the K-dimensional eigenvector space
  phi_g <- drop(crossprod(phi, stat))
  
  ## Fix 2: Add Tikhonov (ridge) regularization to bound the influence of the smallest retained eigenvalues
  Dobs <- sum(phi_g^2 / (lambda_vals + ridge_penalty))
  
  ## Fix 3: Convert natural log to log10 to match your volcano plot axis
  log_p_ln <- pchisq(Dobs, df = K, lower.tail = FALSE, log.p = TRUE) 
  log_p_10 <- log_p_ln / log(10)
  return(-log_p_10)
}
