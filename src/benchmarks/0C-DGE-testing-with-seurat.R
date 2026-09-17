## read the pseudobulk counts deseq object and prepare seurat object with the counts and metadata
## compare results with DESeq differential testing
source('src/scripts/00-load-packages.R')
CONTRAST <- 'old_vs_adult'
# library(DESeq2)
# domain_dds <- read_rds('data/temp/domain_dds.rds')
# domain_dds

## prepare seurat object
# brain <- CreateSeuratObject(counts = counts(domain_dds, normalized = F),
#                             assay = 'RNA', 
#                             meta.data = data.frame(colData(domain_dds)))

#------------------------------------------------------------------------------------------------
## prepare seurat object with domain-level pseudo-bulk counts
#------------------------------------------------------------------------------------------------
library(Seurat)
brain <- read_rds('output/brain_tissue_domains_seurat.rds')
brain

## process the data
brain <- NormalizeData(brain)
brain <- ScaleData(brain)
brain <- RunPCA(brain, features = rownames(brain))

ElbowPlot(brain, ndims = 50)
VizDimLoadings(brain, dims = 1:2)
DimPlot(brain, group.by = c('tissue', 'donor_id', 'age'))

brain <- FindNeighbors(brain, dims = 1:30)
brain <- RunUMAP(brain, dims = 1:30)
DimPlot(brain, group.by = 'age')
genes <- rownames(brain)
write_rds(brain, 'output/brain_tissue_domains_seurat.rds')
#-----------------------------------------------------------------
## perform differential testing with MAST
#-----------------------------------------------------------------
## perform DGE testing
source('src/R/get_cct_neglogp.R')
source('src/R/summarize_dge_results.R')
source('src/R/get_ranked_dge_plot.R')
## iterate over each tissue domain and set the contrast to extract values
mast_dge_res <- list()
for(domain in levels(brain$tissue)){
  message(domain)
  mast_dge_res[[domain]] <- FindMarkers(brain, test.use = 'MAST', latent.vars = 'nCount_RNA',
                                        logfc.threshold = 0, min.pct = 0,
                                        ident.1 = glue('90wk_{domain}'), ident.2 = glue('24wk_{domain}'), 
                                        group.by = 'age_domain')
}

## combine the pvals across domains for each gene
sapply(mast_dge_res, nrow)

## order the genes and combine p-values
mast_res <- summarize_dge_results(mast_dge_res, genes)

p1 <- get_ranked_dge_plot(mast_res$dge)
p2 <- get_ranked_dge_plot(mast_res$dge, statistic = 'stats_min')

pdf(glue('figs/MAST_domain_{CONTRAST}_sorted_stats.pdf'), width = 10, height = 6)
wrap_plots(p1+p2)
dev.off()

write_rds(mast_res, 'output/mast_dge_results.rds')
#-----------------------------------------------------------------
## perform differential testing with DESeq
#-----------------------------------------------------------------
deseq_dge_res <- list()
for(domain in levels(brain$tissue)){
  message(domain)
  deseq_dge_res[[domain]] <- FindMarkers(brain, test.use = 'DESeq2', 
                                        logfc.threshold = 0, min.pct = 0,
                                        ident.1 = glue('90wk_{domain}'), ident.2 = glue('24wk_{domain}'), 
                                        group.by = 'age_domain')
}

## combine the pvals across domains for each gene
sapply(deseq_dge_res, nrow)
sapply(deseq_dge_res, function(res) sum(res['p_val_adj'] < .05, na.rm=T))

## order the genes and combine p-values
deseq_res <- summarize_dge_results(deseq_dge_res, genes)

p1 <- get_ranked_dge_plot(deseq_res$dge, method = 'DESeq2')
p2 <- get_ranked_dge_plot(deseq_res$dge, statistic = 'stats_min', method = 'DESeq2')

pdf(glue('figs/deseq_domain_{CONTRAST}_sorted_stats.pdf'), width = 10, height = 6)
wrap_plots(p1+p2)
dev.off()

write_rds(deseq_res, 'output/deseq_dge_results.rds')


#-----------------------------------------------------------------
## compare the seruat deseq results with direct testing results
#-----------------------------------------------------------------
dx1 <- deseq_res$dge %>% arrange(desc(stats_combined))
dx2 <- vroom(glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv')) %>% 
  mutate(stats_min =  neglogp_min * sign(maxlogFC),
        stats_combined = neglogp_combined * sign(maxlogFC)) %>% 
  arrange(match(gene, dx1$gene)) 

plot(dx1$stats_combined, dx2$stats_combined)
plot(dx1$neglogp_combined~ dx2$neglogp_combined)

#-----------------------------------------------------------------
#-----------------------------------------------------------------


