cct_p.vals <- function(pval, w = NULL){
  if(is.null(w)) 
    w <- rep(1/length(pval), length(pval))
  w <- w / sum(w)
  t_stat <- sum(w * tan((0.5 - pval) * pi))
  if (t_stat > 1e15) {
    log_pcc <- -log(pi) - log(t_stat)
  } else {
    # Standard formula is completely safe for negative, zero, or small positive t_stat
    pcc <- 0.5 - atan(t_stat) / pi
    log_pcc <- log(pcc)
  }
  return(-log_pcc)
}
