library(DESeq2)
library(vsn)

CONTRAST <- 'old_vs_adult'
#--------------------------------
## run pseudo-bulk testing at domain level
#--------------------------------
domain_dds <- read_rds('data/temp/domain_dds.rds')
colData(domain_dds)

domain_dds <- DESeq(domain_dds)
res <- results(domain_dds)
summary(res, alpha = .05)

## gather results across all bins
res.list <- list()
for(domain in levels(domain_dds$tissue)){
  res.list[[domain]] <- results(domain_dds, contrast= c('age_domain', paste('90wk', domain, sep='_'),
                                                  paste('24wk', domain, sep='_')))
}

## get a summary table of p.vals and logFC estimates for all genes across comparisons
res.list <- lapply(res.list, as.data.frame)
write_rds(res.list, glue('data/temp/{CONTRAST}_domain_level_deseq_res.rds'))

#--------------------------------
## gather results across bins
#--------------------------------
source('src/R/get_cct_neglogp.R')

## extract statistics for each contrast
p.vals <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(pvalue)))
logFC <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(log2FoldChange)))
rownames(logFC) <- rownames(p.vals) <- rownames(res.list[[1]])
nlpcx <- apply(p.vals, 1, cct_p.vals)

save(p.vals, logFC, nlpcx, file = glue('output/DESeq_domain_{CONTRAST}_results.rdata'))
#--------------------------------
## compare results with min p-val across domains
#--------------------------------
load('output/DESeq_domain_old_vs_adult_results.rdata')
#--------------------------------
# extract FC of most significant change across bins
#--------------------------------
max.idx <- apply(p.vals, 1, which.min)
max.delta <- logFC[cbind(1:nrow(logFC), max.idx)]
names(max.delta) <- rownames(logFC)
min.neglogp <- -log(p.vals[cbind(1:nrow(p.vals), max.idx)])

## make separate summary dataframe
deseq_dge <- data.frame(neglogp_combined = nlpcx, 
                        maxlogFC = max.delta,
                        neglogp_min = min.neglogp,
                        index = max.idx) %>% 
  mutate(stats =  neglogp_min * sign(maxlogFC)) %>% 
  arrange(desc(stats)) %>% 
  rownames_to_column('gene') %>% 
  mutate(Rank = row_number()) %>%
  as_tibble()

## no. of significant genes
deseq_dge %>%
  dplyr::filter(`neglogp_min` > -log(0.05)) %>%
  with(table(maxlogFC > 0)) ## 102 up/ 94 down

with(deseq_dge, plot(neglogp_combined~maxlogFC, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
with(deseq_dge, plot(neglogp_min~maxlogFC, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
write_csv(deseq_dge, glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv'))


pdf(glue('figs/DESeq_domain_{CONTRAST}_sorted_stats_minlogp.pdf'), width = 7, height = 7)
ggplot(deseq_dge, aes(x = Rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1, alpha = .5) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
             color = "red", size = 2.5, alpha=0.25) +
  
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
    title = glue("Sorted Statistics {CONTRAST} \n domain-level DESeq min-pval")
  ) +
  # coord_cartesian(ylim = c(NA, 20), clip = "off") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12)
  )
dev.off()
