compute_analytical_pvalue <- function(sigma_statistic, test_statistic){
  
  eigen_decomp <- eigen(sigma_statistic, symmetric = TRUE)
  evals <- eigen_decomp$values  
  
  ## restrict to top-K eigenvals
  threshold <- 1e-6*evals[1]
  K <- sum(evals >= threshold)
  
  lambda_vals <- evals[1:K]
  phi <- eigen_decomp$vectors[, 1:K, drop = FALSE]
  
  ## project test statistic vector into the K-dimensional eigenvector space
  phi_g <- drop(crossprod(phi, test_statistic))
  
  Dobs <- sum(phi_g^2 / lambda_vals)
  
  log_p <- pchisq(Dobs, df = K, lower.tail = FALSE, log.p = T) 
  return(-log_p)
}
