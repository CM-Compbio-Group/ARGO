#' **0A-generate_binned_counts** : use multi-gaston output to generate pseudo-bulk counts along isodepth axis with binned data
#' 
#' input : anndata object with multi-gaston inferred isodepth values and raw counts; number of bins to split data into
#' output: raw pseudo-bulk counts in each bin : binned_counts
#'         bin-level metadata : bin_meta with age --> donor --> slide --> bin information, coordinate axes of bins
#' Appends the vst-transformed counts to the bin-metadata, adds reflective padding and saves object
#' 
#' pseudo-bulk counts with either binned or tissue-level resolution is used as input for all differential testing methods
#' continuous DGE tools are input bin-level data from each slide and associated coordinates: funGP, tradeseq and lamian
#' bulk DGE testing tools are input tissue-level data from each slide : deseq2, edgeR and limma-voom



# load libraries
source('src/scripts/00-load-packages.R')
source('src/R/add_reflective_padding.R')
library(anndata)
library(DESeq2)
library(vsn)

# read anndata object
# ad <- read_h5ad('../abundance-modeling/data/raw/integrated_anndata.h5ad')
# obs <- read_rds('data/temp/ad_obs.rds')
ad <- read_h5ad('~/Desktop/integrated_anndata.h5ad')

# bin data based on (binned)-isodepth values
# these bins serve as observed data
# original binned-isodepth coordinates serve as test points for fitting and DGE testing
obs <-  ad$obs %>% rownames_to_column('spot_id') %>% as_tibble() %>% mutate(bin = ntile(binned_isodepth, n = 25)) 

## assign coordinates as mean scaled-isodepth in each bin
obs <- left_join(obs, obs %>% group_by(bin) %>% 
                   summarise(coord_x = mean(binned_isodepth, na.rm=T)),
                 by = 'bin')

x_observed <- obs %>% pull(coord_x) %>% unique() %>% sort()

#--------------------------------
## visualize bins
#--------------------------------
pheatmap::pheatmap(table(obs$bin, obs$tissue), cluster_rows = F)

# distribution of isodepth values across tissue domains by slide
ggplot(obs, aes(x = factor(tissue), y = binned_isodepth, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + # Custom color palette
  theme_classic() +       
  labs(x = "Tissue domain", y = "Isodepth") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) + 
  facet_wrap(~ donor_id)

#----------------------------------------------------------------------
## extract counts data and convert to pseudobulk bin-level data
#----------------------------------------------------------------------
## extract raw counts data 
cts <- ad$layers['counts']
all(rownames(cts) == obs$spot_id)

# define isodepth bin-level counts and metadata
slide_bin_factor <- paste(obs$batch, obs$bin, sep = '_')
binned_counts <- rowsum(as.matrix(cts), group = slide_bin_factor)

## append coordinates of bins
bin_meta <- obs %>% dplyr::count(age, donor_id, batch, coord_x, bin, name = 'n_cells') %>% 
  filter(n_cells > 25) %>%
  mutate(slide_bin = paste(batch, bin, sep = '_'), 
         age_bin = paste(age, bin, sep = '_')) %>%
  mutate(slide_bin = factor(slide_bin, 
                            levels = gtools::mixedsort(unique(slide_bin))),
         age_bin = factor(age_bin, 
                          levels = gtools::mixedsort(unique(age_bin)))
  ) %>% 
  distinct() %>% 
  arrange(slide_bin)

## have some missing bins in slide-12
image(table(bin_meta$bin, bin_meta$batch))
#--------------------------------
## save domain coordinates for downstream functions
x_observed <- sort(unique(bin_meta$coord_x)); (N <- length(x_observed))
x_test <- sort(unique(obs$binned_isodepth)); (M <- length(x_test))
save(x_test, M, file = 'data/temp/test_domain.rdata')
#--------------------------------
##!!!!! ordering by factor orders by the levels !!!!!!##
binned_counts <- binned_counts[as.character(bin_meta$slide_bin), ] 
all(rownames(binned_counts) == bin_meta$slide_bin)
boxplot(binned_counts[, 'C4b']~ bin_meta$age_bin)
write_rds(bin_meta, 'data/temp/bin_meta.rds')
#--------------------------------
## prepare deseq objects to normalize data
#--------------------------------
bin_dds <- DESeqDataSetFromMatrix(countData = t(binned_counts),
                                  colData = bin_meta,
                                  design = ~ age_bin)

bin_dds <- estimateSizeFactors(bin_dds)
hist(sizeFactors(bin_dds))


par(mfrow = c(2,1),  mar = c(3, 3, 4, 2) + 0.1, family = 'Arial', ps=12)
plotCounts(bin_dds, 'C4b', intgroup = 'age_bin')
plotCounts(bin_dds, 'Ctss', intgroup = 'age_bin')
write_rds(bin_dds, 'data/temp/bin_dds.rds')
#--------------------------------
# normalize the counts data
vsd <- varianceStabilizingTransformation(bin_dds, blind = F)
meanSdPlot(assay(vsd))
## save necessary variables (normalized counts and metadata) to single object
save(vsd, bin_meta, file = 'data/temp/binned_data.rdata')
##============================================
## create object with bin-level normalized counts 
## per gene and the metadata information
## add reflective padding to the counts data to assist with fitting GPs
##============================================
# extract normalized counts, and center them before appending to metadata object
genes <- rownames(vsd) # variable names
## append counts data to the metadata object
all(colnames(vsd) == bin_meta$slide_bin)
bin_meta <- bin_meta %>% mutate(bin = as.character(bin)) %>%
  bind_cols(scale(t(assay(vsd)), center = T, scale = F))

## add padding to data
bin_meta <- add_reflective_padding(bin_meta, n_bins = 10, last_bin = 25)
write_rds(bin_meta, 'data/processed/bin_meta.rds')
#--------------------------------------
## prepare tissue-level deseq object
#--------------------------------------
## rename the domains
levels(obs$tissue) <- str_replace_all(levels(obs$tissue),  c(" " = '_',  "/" = '_'))

## prepare pseudo-bulk counts
slide_domain_factor <- paste(obs$batch, obs$tissue, sep = '_')
domain_counts <- rowsum(as.matrix(cts), group = slide_domain_factor)

## don't need the x-coords in this object
domain_meta <- obs %>% dplyr::count(age, donor_id, batch, tissue, name = 'n_cells') %>% 
  filter(n_cells > 25) %>%
  mutate(slide_domain = paste(batch, tissue, sep = '_'), 
         age_domain = paste(age, tissue, sep = '_')) %>%
  mutate(slide_domain = factor(slide_domain, 
                               levels = gtools::mixedsort(unique(slide_domain))),
         age_domain = factor(age_domain, 
                             levels = gtools::mixedsort(unique(age_domain)))
  ) %>% 
  distinct() %>% 
  arrange(slide_domain)

domain_meta

## plenty of slides don't have ventricle regions
table(domain_meta$tissue, domain_meta$batch)
##!!!!! ordering by factor orders by the levels !!!!!!##
domain_counts <- domain_counts[as.character(domain_meta$slide_domain), ] 
all(rownames(domain_counts) == domain_meta$slide_domain)
boxplot(domain_counts[, 'C4b']~ domain_meta$age_domain)
write_rds(domain_meta, 'data/temp/domain_meta.rds')
#--------------------------------
## prepare the deseq object
domain_dds <- DESeqDataSetFromMatrix(countData = t(domain_counts),
                                     colData = domain_meta,
                                     design = ~ age_domain)
domain_dds <- estimateSizeFactors(domain_dds)
hist(sizeFactors(domain_dds))
plotCounts(domain_dds, 'C4b', intgroup = 'age_domain')
write_rds(domain_dds, 'data/temp/domain_dds.rds')