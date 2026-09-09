get_ranked_dge_plot <- function(df, statistic = 'neglogp',CONTRAST = 'old_vs_adult', method = 'MAST', markers){
  
  df <- df %>% arrange(desc(!!sym(statistic))) %>% mutate(Rank = row_number())
  p <- ggplot(df, aes(x = Rank, y = !!sym(statistic))) +
    # Draws the main sorted line/points
    geom_line(color = "grey70") +
    geom_point(color = "black", size = 1, alpha = .5) +
    
    # Highlight the specific genes with a unique color
    # geom_point(data = . %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
    #            color = "red", size = 2.5, alpha=0.25) +
    
    geom_point(data = . %>% filter(gene %in% markers), color = 'red',size=2.5, alpha=0.25 ) + 
    # Add the non-overlapping labels
    geom_text_repel(
      # data = .  %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
      data = . %>% filter(gene %in% markers),
      aes(label = gene),
      size = 4,
      color = "black",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "red", 
      max.overlaps = 50
    ) +
    
    # Labels and clean styling
    labs(
      x = "Rank",
      y = "Statistic Value",
      title = glue("Sorted Statistics {CONTRAST} \n {method}")
    ) +
    # coord_cartesian(ylim = c(NA, 20), clip = "off") + 
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12)
    )
  return(p)
}

