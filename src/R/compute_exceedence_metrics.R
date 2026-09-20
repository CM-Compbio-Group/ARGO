compute_exceedence_metrics <- function(res, x_test){
  stat <- res[['statistic']]
  band <- res[['band']]
  M <- length(x_test)
  
  
  ## find max deviation outside the conf bands
  max_dev <- ifelse(any(abs(stat) > band), max(abs(stat)), 0)
  ## compute directional stats 
  poe_up <- mean(stat > band)
  poe_dn <- mean(stat < -band)
  
  dx <- diff(x_test)
  ## compute area relative to zero
  exc_up <- (pmax(stat - band, 0) > 0)
  stat_exc_up <- ifelse(exc_up > 0 ,stat, 0)
  
  auec_up <- sum((stat_exc_up[-1] + stat_exc_up[-M]) / 2 * dx) # auc with mean b/w adjacent points
  
  exc_dn <- pmax(-stat - band, 0) #flipping the curve
  stat_exc_dn <- ifelse(exc_dn > 0, -stat, 0) # compute absolute area in negative direction
  auec_dn <- sum((stat_exc_dn[-1] + stat_exc_dn[-M]) / 2 * dx)
  
  return(c(
     'max_deviation' = max_dev,
     'poe_up' = poe_up,
     'poe_down' = poe_dn,
     'auec_up' = auec_up,
     'auec_down' = auec_dn))
}


