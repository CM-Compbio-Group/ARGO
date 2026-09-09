generate_cc_matrix <- function(contrast_vector){
  kk <- length(contrast_vector) ## number of latent vectors
  index <- which(contrast_vector == 1) ## where the one lies
  cc_base <- diag(kk - 1)
  cc_vector_rem <- matrix(-contrast_vector[-index], ncol = 1) 
  
  if(index == 1){
    cc <- cbind(cc_vector_rem, cc_base)
  } else if (index == kk){
    cc <- cbind(cc_base, cc_vector_rem)  
  } else{
    cc <- cbind(cc_base[, 1:(index - 1), drop = FALSE], 
                cc_vector_rem,
                cc_base[, index : (kk - 1), drop = FALSE])
  }
  return(cc)
}