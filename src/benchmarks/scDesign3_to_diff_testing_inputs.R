#' script to go from simulated counts using scDesign3 to inputs for various differential testing methods

## load counts data
load('data/scdesign_simulation.rdata')

## prepare pseudobulk counts and prepare deseq object
# prepare bins for observations and testing. 
# The 2087 values indicates individual cells. don't need to fit to all cells
sim_meta <-  sim_meta %>% 
  mutate(bin_observed = factor(ntile(pseudotime, n = 25)), ## 25 points where data is observed
         bin_test = factor(ntile(pseudotime, n = 100)), # 100 points to test fits at
         sample_bin = paste(sample, bin_observed, sep = '_'),
         condition_bin = paste(condition, bin_observed, sep = '_'))

sim_meta <- sim_meta %>% mutate(condition_bin = factor(condition_bin, 
                                                       levels =  gtools::mixedsort(unique(condition_bin))),
                                sample_bin = factor(sample_bin, 
                                                    levels = gtools::mixedsort(unique(sample_bin)))) 
## append coordinates for each bin
sim_meta <- right_join(sim_meta, sim_meta %>% group_by(bin_observed) %>% 
                         summarise(coord_x = mean(pseudotime)), by = 'bin_observed')
sim_meta <- right_join(sim_meta, sim_meta %>% group_by(bin_test) %>% 
                         summarise(test_x = mean(pseudotime)), by = 'bin_test')
sim_meta

write_rds(sim_meta, file = 'output/simulations/sim_meta.rds')

## prepare bin level metadata
bin_meta <- sim_meta %>% dplyr::count(condition, replicate, sample, 
                                      bin_observed, coord_x, condition_bin, sample_bin,
                                      name = 'n_cells') %>% 
  distinct() %>% 
  arrange(sample_bin)

## convert counts to pseudobulk values
binned_counts <- rowsum(t(sim_counts), group = sim_meta$sample_bin)
all(rownames(binned_counts) == bin_meta$sample_bin)
#===================================================================
## prepare deseq object with pooled bin-level counts data
#===================================================================
dds <- DESeqDataSetFromMatrix(countData = t(binned_counts),
                              colData = bin_meta,
                              design = ~ 0 + condition_bin)
dds <- estimateSizeFactors(dds)

## visualize counts for some genes
plotCounts(dds, gene = 'Npepl1', intgroup = 'bin_observed')
plotCounts(dds, gene = 'Ckb', intgroup = 'bin_observed')

## perform counts normalization
vsd <- varianceStabilizingTransformation(dds, blind = F)
plotPCA(vsd, intgroup = 'bin_observed') 

save(dds, vsd, bin_meta, file = 'output/simulations/deseq_objects.rdata')


#===================================================================
## append the counts to bin meta to prepare input for funGP
#===================================================================
## append the counts data to the bin_meta object
all(colnames(vsd) == bin_meta$sample_bin)
## append centered log-normalized counts data
bin_meta <- bin_meta %>%
  # mutate(bin = as.character(bin)) %>%
  bind_cols(scale(t(assay(vsd)), center = T, scale = F))

plot(bin_meta$Chgb~ bin_meta$condition_bin)
plotCounts(dds, gene = 'Chgb', intgroup = 'condition_bin')

## get variables for fitting continuous-time models
x_observed <- sort(unique(bin_meta$coord_x)); N <- length(x_observed)
x_test <- sort(unique(sim_meta$test_x)); M <- length(x_test)
save(bin_meta, x_observed, x_test, N, M, file = 'output/simulations/fungp_inputs.rdata')
