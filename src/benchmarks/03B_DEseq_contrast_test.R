library(DESeq2)
library(vsn)

CONTRAST <- 'old_vs_adult'
#--------------------------------
## run pseudo-bulk testing at bin level
#--------------------------------
bin_dds <- read_rds('data/temp/bin_dds.rds')
colData(bin_dds)

bin_dds <- DESeq(bin_dds)
res <- results(bin_dds)
summary(res, alpha = .05)

## gather results across all bins
res.list <- list()
for(bin in levels(bin_dds$bin)){
  res.list[[bin]] <- results(bin_dds, contrast= c('age_bin', paste('90wk', bin, sep='_'),
                                                  paste('24wk', bin, sep='_')))
}

## get a summary table of p.vals and logFC estimates for all genes across comparisons
res.list <- lapply(res.list, as.data.frame)
write_rds(res.list, glue('data/temp/{CONTRAST}_bin_level_deseq_res.rds'))

#--------------------------------
## gather results across bins
#--------------------------------
source('src/R/get_cct_neglogp.R')

## extract statistics for each contrast
p.vals <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(pvalue)))
logFC <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(log2FoldChange)))
rownames(logFC) <- rownames(p.vals) <- rownames(res.list[[1]])
nlpcx <- apply(p.vals, 1, cct_p.vals)

save(p.vals, logFC, nlpcx, file = glue('output/DESeq_{CONTRAST}_results.rdata'))
# extract FC of most significant change across bins
max.idx <- apply(p.vals, 1, which.min)
max.delta <- logFC[cbind(1:nrow(logFC), max.idx)]
names(max.delta) <- rownames(logFC)

## make separate summary dataframe
deseq_dge <- data.frame(neglogp_combined = nlpcx, 
                       maxlogFC = max.delta,
                       index = max.idx) %>% 
  mutate(stats =  neglogp_combined * sign(maxlogFC)) %>% 
  arrange(desc(stats)) %>% 
  rownames_to_column('gene') %>% 
  mutate(Rank = row_number()) %>%
  as_tibble()

## no. of significant genes
deseq_dge %>%
  dplyr::filter(`neglogp_combined` > -log(0.05)) %>%
  with(table(maxlogFC > 0)) ## 86 up/ 86 down

with(deseq_dge, plot(neglogp_combined~maxlogFC, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
write_csv(deseq_dge, glue('data/sheets/DESeq_DGE_bin_{CONTRAST}.csv'))


pdf(glue('figs/DESeq_bin_{CONTRAST}_sorted_stats.pdf'), width = 8, height = 8)
ggplot(deseq_dge, aes(x = Rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
             color = "red", size = 2.5, alpha=0.5) +
  
  # Add the non-overlapping labels
  geom_text_repel(
    data = .  %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
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
    title = glue("Sorted Statistics {CONTRAST} \n bin-level DESeq combined-pvals")
  ) +
  # coord_cartesian(ylim = c(NA, 20), clip = "off") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12)
  )
dev.off()
