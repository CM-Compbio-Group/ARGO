## KLT approximation of posterior covariance
KLT_CI_bands <- function(sigma, n_samples = 1e3, conf.level = .99){
  
  ## eigen decomposition of the cov kernel under H0
  eig <- eigen(sigma)
  ## filter by eigen values
  lmax <- eig$values[1]
  m <- which(eig$values < 1e-6*lmax)[1]
  
  if(is.na(m)){
    m <- length(eig$values)
  }
  
  ## represent the process as sum of iid std normal variables
  lambda <- diag(eig$values[1:m])
  U <- eig$vectors[, 1:m]
  
  ## not worrying about the non-central chi-squared parameter
  threshold <- qchisq(conf.level, m) # quantile for chi-sq distribution w/ 'm' DOF
  
  ## no. of samples to draw with expected acceptance rate
  n_to_draw <- ceiling(n_samples / conf.level * 1.25)
  
  Z_raw <- matrix(rnorm(n_to_draw * m), nrow = m)
  
  ## which samples lie within conf level
  sq_norms <- colSums(Z_raw^2)
  Z_accepted <- Z_raw[, sq_norms < threshold, drop = FALSE]
  
  Z <- Z_accepted[, 1:n_samples]
  
  ## compute the process
  G = U %*% sqrt(lambda) %*% Z
  
  ## find the limits of the confidence band using quantile estimate for more stability
  band <- apply(G, 1, quantile, .999)
  
  return(list('band' = band, 'dof' = m))  
  
}
