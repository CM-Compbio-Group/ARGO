get_pc_plot <- function(vsd, intgroup){
  # 1. Extract the PCA data rather than drawing the plot directly
  pcaData <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  
  # 2. Calculate the variance percentages for the axis labels
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # 3. Build the custom plot with your desired dot size (e.g., size = 1.5)
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!sym(intgroup))) +
    geom_point(size = 1.5, alpha = .5) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(color = intgroup) +
    theme_bw() # Optional clean theme
  
  return(p)
}