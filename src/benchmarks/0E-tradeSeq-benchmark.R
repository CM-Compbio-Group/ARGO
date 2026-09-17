BiocManager::install("tradeSeq")
library(tradeSeq)


## load bin level counts data and metadata
bin_dds <- read_rds('data/temp/bin_dds.rds')
bin_meta <- read_rds('data/temp/bin_meta.rds')
all(colnames(bin_dds) == bin_meta$slide_bin)

## select no. of knots for GAM fitting
## use the isodepth values for each bin as pseudotime

icMat <- evaluateK(counts = counts(bin_dds, normalized = F),
                   pseudotime = bin_meta$coord_x, 
                   cellWeights = rep(1, nrow(bin_meta)),
                   conditions = bin_meta$age, 
                   nGenes = 374,
                   k = 3:10)


## fit GAM and get an sce object
sce <- fitGAM(counts = counts(bin_dds, normalized = F), 
              pseudotime = bin_meta$coord_x, 
              cellWeights = rep(1, nrow(bin_meta)),
              conditions = bin_meta$age, 
              nknots = 7, 
              verbose = T, 
              parallel = T)

mean(rowData(sce)$tradeSeq$converged)
sce


plotSmoothers(sce, assays(sce)$counts, gene = "C4b", alpha = 1, border = TRUE) + ggtitle("C4b")
plotSmoothers(sce, assays(sce)$counts, gene = "Ctss", alpha = 1, border = TRUE) + ggtitle("Ctss")
plotSmoothers(sce, assays(sce)$counts, gene = "Grin2b", alpha = 1, border = TRUE) + ggtitle("Grin2b")
plotSmoothers(sce, assays(sce)$counts, gene = "Ctsd", alpha = 1, border = TRUE) + ggtitle("Ctsd")


## run DGE test between conditions
condRes <- conditionTest(sce, pairwise = T)
mean(condRes$pvalue_conds2vs3 <= 0.05, na.rm = TRUE)

## order genes by significance
condRes <- condRes %>% rownames_to_column('gene') %>% arrange(desc(waldStat_conds2vs3))
plot(seq(374), condRes$waldStat_conds2vs3)

oo <- order(condRes$waldStat_conds2vs3, decreasing = TRUE)

plotSmoothers(sce, assays(sce)$counts,
              gene = rownames(assays(sce)$counts)[oo[106]],
              alpha = 1, border = TRUE)


write_rds(condRes, file = 'output/tradeSeq_res.rds')
