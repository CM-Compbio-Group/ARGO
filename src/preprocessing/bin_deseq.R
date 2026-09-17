source('src/R/get_pc_plot.R')
p1 <- get_pc_plot(vsd, 'donor_id')
p2 <- get_pc_plot(vsd, 'bin')
p3 <- get_pc_plot(vsd, 'age')

combined_plot <- wrap_plots(list(p3, p1, p2), ncol=3) + 
  plot_layout(guides = "keep")

ggsave(combined_plot, 
       filename = file.path('figs', glue('bin_counts_PCA.pdf')),
       width=12,height=4)

#--------------------------------
## run pseudo-bulk testing at bin level
#--------------------------------
bin_dds <- DESeq(bin_dds)
res <- results(bin_dds)
summary(res, alpha = .05)

## gather results across all bins
res.list <- list()
for(bin in levels(bin_dds$bin)){
  res.list[[bin]] <- results(bin_dds, contrast= c('age_bin', paste('90wk', bin, sep='_'),
                                                  paste('4wk', bin, sep='_')))
}

## get a summary table of p.vals and logFC estimates for all genes across comparisons
res.list <- lapply(res.list, as.data.frame)
write_rds(res.list, 'data/temp/90wk_vs_4wk_bin_level_deseq_res.rds')
#--------------------------------
## get cumulative results from the test
#--------------------------------
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

## extract statistics for each contrast
p.vals <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(pvalue)))
logFC <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(log2FoldChange)))
rownames(logFC) <- rownames(p.vals) <- rownames(bin_dds)
nlpcx <- apply(p.vals, 1, cct_p.vals)

# extract FC of most significant change across bins
max.idx <- apply(p.vals, 1, which.min)
max.delta <- logFC[cbind(1:nrow(logFC), max.idx)]
names(max.delta) <- rownames(logFC)

## make separate summary dataframe
p.val.df <- data.frame(neglogp_combined = nlpcx, 
                       maxlogFC = max.delta,
                       index = max.idx) %>% 
  mutate(stats =  neglogp_combined * sign(maxlogFC)) %>% 
  arrange(desc(stats)) %>% 
  rownames_to_column('gene') %>% 
  mutate(rank = row_number()) %>%
  as_tibble()

## no. of significant genes
p.val.df %>%
  dplyr::filter(`neglogp_combined` > -log(0.05)) %>%
  with(table(maxlogFC > 0)) ## 86 up/ 86 down

with(p.val.df, plot(neglogp_combined~maxlogFC, pch=20, col = adjustcolor("red", alpha.f = 0.5)))

write_csv(p.val.df, 'data/sheets/DESeq_DGE_bin_old_vs_juv.csv')
#--------------------------------------------
## visualize results for the top genes
#--------------------------------------------
to_plot <- p.val.df %>% pull(gene) %>% {c(head(., 20), tail(., 20))}
plotCounts(bin_dds, 'C4b', intgroup = 'age_bin', pch=20)

pdf('figs/DEseq_bin_LFC_topDiff_heatmap.pdf', width = 8, height = 8)
Heatmap(logFC[to_plot, ], name = 'LFC', 
        show_row_names = T,
        cluster_columns = F, cluster_rows = F,
        column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
        heatmap_legend_param = list(
          title = "LogFC",
          title_gp = gpar(fontsize = 16, fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
          legend_height = unit(4, "cm"), # Adjust height of the color bar
          grid_width = unit(5, "mm")      # Adjust width of the color bar
        ))
dev.off()

#--------------------------------------------
## make the ordered stats plot
#--------------------------------------------
pdf('figs/DESeq_bin_sorted_stats.pdf', width = 8, height = 8)
ggplot(p.val.df, aes(x = rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% filter(gene %in% to_plot),
             color = "red", size = 2.5, alpha=0.5) +
  
  # Add the non-overlapping labels
  geom_text_repel(
    data = . %>% filter(gene %in% to_plot),
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
    title = "Sorted Statistics Profile with \n bin-level DESeq combined-pvals"
  ) +
  # coord_cartesian(ylim = c(NA, 20), clip = "off") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12)
  )
dev.off()
