## compute statistics from funGP curves and conf bands

compute_stat_binary <- function(exact_pvals, M){
  stat_binary <- lapply(exact_pvals, function(res) {
    stat <- res[['statistic']]; 
    band <- res[['band']];
    case_when(
      stat < -band ~ -1,
      stat > band ~ 1,
      TRUE ~ 0
    )
  })
  
  stat_binary <- do.call('rbind', stat_binary)
  colnames(stat_binary) <- paste0('bin', 1:M)
  return(stat_binary)
}
  
compute_stat_continuous <- function(exact_pvals, M){
  stat_lfc <- lapply(exact_pvals, function(res) {
    stat <- res[['statistic']]; 
    band <- res[['band']];
    case_when(
      abs(stat) > band ~ (stat - band),
      TRUE ~ 0
    )
  })
  
  stat_lfc <- do.call('rbind', stat_lfc)
  colnames(stat_lfc) <- paste0('bin', 1:M)
  return(stat_lfc)
}
