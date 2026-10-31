summarize_argo_excursions <- function(obj){
  s <- obj$statistic
  b <- obj$band

  
  pos_roi <- which(s > b)  
  neg_roi <- which(s < -b)
  
  pos_mass <- sum(s[pos_roi])
  neg_mass <- sum(abs(s[neg_roi]))
  
  net_deviation <- pos_mass - neg_mass
  net_sign <- ifelse(net_deviation > 0, 1, -1)
  polarity <- max(pos_mass, neg_mass)/pmax(pos_mass + neg_mass, 1e-3)
  tibble(neglogp = obj$neglogp, 
         raw_max_dev = obj$deviation, 
         max_pos_dev = ifelse(any(pos_roi), max(s[pos_roi]), max(s)),
         max_neg_dev = ifelse(any(neg_roi), min(s[neg_roi]), min(s)),
         pos_mass = pos_mass,
         neg_mass = neg_mass,
         n_pos_bins = length(pos_roi),
         n_neg_bins = length(neg_roi),
         net_sign = net_sign, 
         directed_deviation = ifelse(net_sign > 0, max_pos_dev, max_neg_dev),
         polarity = polarity
         )
}
