bin_kdtree_adaptive <- function(data, target_spot_count = 100) {
  # data: matrix or data.frame with 4 columns (spots as rows)
  # target_spot_count: target number of spots per bin
  
  n_obs <- nrow(data)
  
  # Helper function for recursive conditional median splitting
  split_node <- function(indices, depth = 0) {
    n_node <- length(indices)
    
    # Base case: stop splitting when target leaf size is reached
    if (n_node <= target_spot_count * 1.5) {
      return(rep(1, n_node))
    }
    
    # Pick split axis: choose dimension with maximum variance among local points
    sub_data <- data[indices, , drop = FALSE]
    vars <- apply(sub_data, 2, var)
    
    # Fallback to cycling axes if variance is zero/flat
    split_dim <- if (all(vars == 0)) (depth %% ncol(data)) + 1 else which.max(vars)
    
    # Order local indices along the chosen axis to guarantee exact 50/50 split
    ord <- indices[order(sub_data[, split_dim])]
    half <- floor(n_node / 2)
    
    left_idx  <- ord[1:half]
    right_idx <- ord[(half + 1):n_node]
    
    left_bins  <- split_node(left_idx, depth + 1)
    right_bins <- split_node(right_idx, depth + 1)
    
    # Return mapping
    out <- numeric(n_node)
    out[match(left_idx, indices)]  <- left_bins
    out[match(right_idx, indices)] <- right_bins + max(left_bins)
    return(out)
  }
  
  bin_assignments <- split_node(1:n_obs)
  return(factor(bin_assignments))
}