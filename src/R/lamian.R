## run differential testing for the {CONTRAST} using Lamian

library(Lamian)
library(DESeq2)

## load the normalized counts data
## lamian runs with log-transformed counts as input
bin_dds <- read_rds('data/temp/bin_dds.rds')
load('data/temp/input_binned_data.rdata') # loads the variance-stabilized counts and associated metadata
all(bin_meta$slide_bin == colnames(vsd))

## pseudotime values
pt <- bin_meta$coord_x; names(pt) <- bin_meta$slide_bin

# set 24wk as baseline age to contrast with 90wk
sample_meta <- bin_meta %>% select(age, batch) %>% 
  distinct() %>% mutate(age = fct_relevel(age, '24wk'))
sample_meta

design <- model.matrix(~age, data = sample_meta)
rownames(design) <- sample_meta$batch
design


# get results for old-vs-adult contrast controlling for donor effects
res <- lamian_test(expr = as.matrix(assay(vsd)),
                   cellanno = bin_meta %>% select(Cell = slide_bin, 
                                                  Sample = batch) %>%
                     as.data.frame(),
                   pseudotime = pt,
                   design = design, 
                   test.type= 'variable', # testing for diff expression with covariate
                   testvar = 3, # test 90wk vs 24wk contrast
                   permuiter =100, 
                   overall.only = T,
                   ncores = 8) 

## get differential dynamic genes statistics
stat <- res$statistics %>% data.frame() %>% 
  mutate(neglogp = -log.pval.overall) %>% 
  arrange(desc(neglogp))

## visualize volcano plot
with(stat, plot(neglogp~ z.overall, col = adjustcolor('red', alpha = .5)))
(diff_genes <- filter(stat, pval.overall < .05) %>% rownames())

plotCounts(bin_dds, 'Ctss', intgroup = 'batch')

#====================================================#
## save results
write_rds(res, 'output/lamian_DGE_res.rds')
