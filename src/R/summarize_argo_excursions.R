summarize_argo_excursions <- function(obj){
  s <- obj$statistic
  b <- obj$band
  
  pos_roi <- which(s > b)  
  neg_roi <- which(s < -b)
  
  raw_max_pos <- max(s)
  raw_max_neg <- min(s)
  
  sig_pos_peak <- ifelse(any(pos_roi), max(s[pos_roi]), 0)
  sig_neg_peak <- ifelse(any(neg_roi), abs(min(s[neg_roi])), 0)
  
  mean_sig_pos <- ifelse(any(pos_roi), mean(s[pos_roi]), 0)
  mean_sig_neg <- ifelse(any(neg_roi), mean(s[neg_roi]), 0)
  
  # net_sign <- ifelse(mean_sig_pos > abs(mean_sig_neg), 1, -1)
  net_sign <- ifelse(sig_pos_peak > sig_neg_peak, 1, -1)
  # net_sign <- ifelse(raw_max_pos > abs(raw_max_neg), 1, -1)
  
  pos_mass <- sum(s[pos_roi])
  neg_mass <- abs(sum(s[neg_roi]))
  
  polarity <- max(pos_mass, neg_mass)/pmax(pos_mass + neg_mass, 1e-8)
  
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

# plot(s, ylim = range(c(s, b, -b)), col = 'blue', pch=20, cex=1.5)
# points(b, col = 'red', pch = 20)
# points(-b, col = 'red', pch = 20)
