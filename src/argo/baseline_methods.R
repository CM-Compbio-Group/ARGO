#' prepare brain aging data for fitting baseline differential testing methods
#' running baseline methods with pseudo-bulk data
#' group at level of cell type annotations -- cell types themselves change with age 
#' group at scale of domain annotations -- more stable representation of data
#' Lamian is run at single cell resolution
#' tradeSeq is run at scale of binned isodepth -- does not scale to fitting at single cell resolution. I'll prepare the data for it separately


ad <- read_h5ad('~/Desktop/integrated_anndata.h5ad')
counts <- as.matrix(ad$layers['counts'])


metadata <-  ad$obs %>% rownames_to_column('cell_id') %>% 
  as_tibble() %>% select(cell_id,  donor_id, sample_id = batch,
                         condition = age, tissue,
                         pseudotime = binned_isodepth) %>% 
  mutate(tissue = str_replace_all(tissue, c(' ' = '_', '-' = '_', '/' = '_'))) %>%
  mutate(condition = factor(condition, levels = c('4wk', '24wk', '90wk')),
         donor_id = factor(donor_id, levels = gtools::mixedsort(unique(donor_id))), 
         tissue = factor(tissue)) %>% 
  mutate(grid_key = paste(sample_id, tissue, sep = '___'), 
         condition_key = paste(condition, tissue, sep = '___')) %>%
  mutate(grid_key = factor(grid_key, 
                           levels = gtools::mixedsort(unique(grid_key))),
         condition_key = factor(condition_key, 
                                levels = gtools::mixedsort(unique(condition_key))))

stopifnot('counts ordering does not match with metadata' = all(rownames(counts) == metadata$cell_id))

# prepare domain level counts data
domain_cts <- rowsum(counts, group = metadata$grid_key)
domain_meta <- metadata %>% select(-c(cell_id, pseudotime)) %>% 
  distinct() %>% arrange(grid_key)
stopifnot('binned counts does not match with metadata' = all(rownames(domain_cts) == domain_meta$grid_key))

#====================================
# filter out bins with high sparsity
#====================================
domain_meta <- domain_meta %>% 
  left_join(dplyr::count(metadata, grid_key), by = 'grid_key') %>%
  dplyr::filter(n > 25)

domain_cts <- domain_cts[domain_meta$grid_key, ]
#====================================
# pass the data for differential testing
#====================================
## prepend X to condition key
domain_meta$condition_key <- paste0('X', domain_meta$condition_key)
diff.testing.methods <- list(
  DESeq2 = run_deseq,
  edgeR = run_edger, 
  Limma = run_limma
)

## call differential testing passing the local environment
res.list <- map(diff.testing.methods, ~ .x(domain_cts, domain_meta))

res.list[['Lamian']] <- run_lamian(counts, metadata)
res.list[['tradeSeq']] <- run_tradeseq(counts, metadata)

res.ARGO <- read_rds('output/old_vs_adult_N30_diff_testing_results.rds')
funGP.df <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())

res.list[['ARGO']] <- funGP.df

write_rds(res.list, file = glue('output/{CONTRAST}_DGE_results_list.rds'))
#===============================================================
## benchmarking
#===============================================================
require(readxl)
require(UpSetR)
require(VennDiagram)
library(PRROC)

## summarize results with upset-plot
res.sig <- lapply(res.list, function(res) filter(res, fdr < 0.05) %>% pull(gene))
sapply(res.sig, length)

pdf('figs/upsetPlot_FDR_05.pdf', width = 8, height = 6)
print(upset(fromList(res.sig), order.by = "freq", nsets = length(res.sig), nintersects = 40), newpage= F)
dev.off()
