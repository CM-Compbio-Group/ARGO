add_reflective_padding <- function(data, n_bins, last_bin = 25){
  #' function to add a reflective padding to observed data
  #' 
  #' @params n_bins : number of bins to extend the boundaries by
  #' @params last_bin : ID of right-boundary in data
  #' 
  #' returns the data with left- and right-paddings appended
  
  coord_x <- unique(data$coord_x)
  delta <- median(diff(coord_x)) 
  ## append left padding
  for(bx in seq(n_bins)){
    bin_y <- filter(data, bin == bx) %>% 
      mutate(bin = glue('L{bx}'),
             coord_x = -delta * bx,
             slide_bin = paste(batch, glue('L{bx}'), sep = '_'),
             age_bin = paste(age, glue('L{bx}'), sep = '_'))
    
    data <- rbind(data, bin_y)
  }
  
  
  ## append right padding
  for(bx in seq(n_bins)){
    bin_y <- filter(data, bin == (last_bin - bx + 1)) %>% 
      mutate(bin = glue('R{bx}'),
             coord_x = 1 + (delta * bx),
             slide_bin = paste(batch, glue('R{bx}'), sep = '_'),
             age_bin = paste(age, glue('R{bx}'), sep = '_'))
    
    data <- rbind(data, bin_y)
  }
  
  data <- data %>% mutate(bin = factor(bin, levels = c(paste0('L', rev(seq(n_bins))),
                                                       1:last_bin,
                                                       paste0('R', seq(n_bins))))) %>%
    arrange(age, donor_id, batch, bin)
  
  return(data)
}