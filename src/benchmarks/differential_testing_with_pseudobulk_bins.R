#' load the scDesign3-simulated counts and perform differential testing with various approaches
#' I need to convert the data into pseudobulk binned counts first, do library size normalization and then feed to the algorithms

library(DESeq2)
library(edgeR)
library(limma)
library(tradeSeq)
library(Lamian)
library(glue)
require(patchwork)

## scripts to gather DGE testing results
source('src/R/get_cct_neglogp.R')
source('src/R/summarize_dge_results.R')
source('src/R/get_ranked_dge_plot.R')

load('output/simulations/deseq_objects.rdata')
#===================================================================
## run differential testing with deseq
#===================================================================
design(dds) <- ~ 0 + condition_bin
dds <- DESeq(dds)
res <- results(dds)
summary(res, alpha = .05)

## gather results across all bins
res.list <- vector("list", 25)
for(bin in seq_len(25)){
  message(bin)
  res.list[[bin]] <- results(dds, contrast= c('condition_bin', paste('condition_B', bin, sep='_'),
                                                        paste('condition_A', bin, sep='_')))
}


head(res.list[[1]])
deseq_res <- summarize_dge_results(res.list, genes = rownames(dds), 
                                   pval = 'pvalue', lfc = 'log2FoldChange')

## summarize findings
p1 <- get_ranked_dge_plot(deseq_res$dge, method = 'DESeq2', markers = effect_genes)
with(deseq_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
deseq_res$dge %>% filter(fdr < 0.05) %>% with(table(max.lfc > 0))
deseq_res$dge %>% filter(fdr < 0.05) %>% with(table(gene %in% null_genes)) ## how many control genes are being called differential
write_rds(deseq_res, 'output/simulations/deseq_DGE_res.rds')
#===================================================================
## run differential testing with edgeR
#===================================================================
y <- DGEList(counts = t(binned_counts), 
             samples = bin_meta)

## normalize the counts
y <- calcNormFactors(y)
## design matrix for contrasts
design <- model.matrix(~ 0 + condition_bin, data = y$samples)
colnames(design) <- gsub("^condition_bin", "", colnames(design))
print(colnames(design))

## run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

## extract results with contrasts
res.list <- list()
for(bin in seq_len(25)){
  message(bin)
  contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
  con <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  ## run the test and extract results
  qlf <- glmQLFTest(fit, contrast = con)
  ## save results to list
  res.list[[bin]] <- topTags(qlf, n = Inf, sort.by = "none")$table
}

names(res.list[[1]])
edger_res <- summarize_dge_results(res.list, genes = rownames(dds), pval = 'PValue', lfc = 'logFC')

p2 <- get_ranked_dge_plot(edger_res$dge, method = 'edgeR', markers = effect_genes)
with(edger_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
edger_res$dge %>% filter(fdr < 0.05) %>% with(table(max.lfc > 0))
edger_res$dge %>% filter(fdr < 0.05) %>% with(table(gene %in% null_genes)) ## how many control genes are being called differential
write_rds(edger_res, file = 'output/simulations/edgeR_DGE_res.rds')

#===================================================================
## run differential testing with limma
#===================================================================
## apply voom-transformation to generate the precision weights
v <- voom(y, design = design, plot = T)
# fit initial linear model
fit <- lmFit(v, design)

## extract results with contrasts
res.list <- vector("list", 25)
for(bin in seq_len(25)){
  message(bin)
  contrast_formula <- glue("condition_B_{bin} - condition_A_{bin}")
  con <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  ## fit contrast to the linear model
  fit_contrast <- contrasts.fit(fit, contrasts =con)
  
  ## run empirical bayes moderation
  fit_ebayes <- eBayes(fit_contrast)
  
  ## save results to list
  res.list[[bin]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
}


limma_res <- summarize_dge_results(res.list, rownames(dds), pval = 'P.Value', lfc = 'logFC')

(p3 <- get_ranked_dge_plot(limma_res$dge, method = 'limma-voom', markers = effect_genes))
with(limma_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
limma_res$dge %>% filter(fdr < .05) %>% with(table(max.lfc > 0)) ## 103 upreg; 87 downreg
limma_res$dge %>% filter(fdr < .05) %>% with(table(gene %in% null_genes)) 
write_rds(limma_res, file = 'output/simulations/limma_DGE_res.rds')

pdf('figs/diff_expression_testing_on_simulations.pdf', width = 15, height=6)
plot(wrap_plots(p1 | p2 | p3))
dev.off()
#===================================================================
## run differential testing with tradeSeq
#===================================================================
## select no. of knots for GAM fitting
## use the isodepth values for each bin as pseudotime
icMat <- evaluateK(counts = counts(dds, normalized = F),
                   pseudotime = bin_meta$coord_x, 
                   cellWeights = rep(1, nrow(bin_meta)),
                   conditions = bin_meta$condition, 
                   nGenes = G,
                   k = 3:10)


## fit GAM and get an sce object
sce <- fitGAM(counts = counts(dds, normalized = F), 
              pseudotime = bin_meta$coord_x, 
              cellWeights = rep(1, nrow(bin_meta)),
              conditions = bin_meta$condition, 
              nknots = 6, 
              verbose = T, 
              parallel = T)

mean(rowData(sce)$tradeSeq$converged)

plotSmoothers(sce, assays(sce)$counts, gene = "Pax6", alpha = 1, border = TRUE) + ggtitle("C4b")
plotSmoothers(sce, assays(sce)$counts, gene = "Cdk4", alpha = 1, border = TRUE) + ggtitle("Ctss")
plotSmoothers(sce, assays(sce)$counts, gene = "Chga", alpha = 1, border = TRUE) + ggtitle("Grin2b")
plotSmoothers(sce, assays(sce)$counts, gene = "Pyy", alpha = 1, border = TRUE) + ggtitle("Ctsd")
plotSmoothers(sce, assays(sce)$counts, gene = "Neurog3", alpha = 1, border = TRUE) + ggtitle("Ctsd")

## run DGE test between conditions
condRes <- conditionTest(sce)
mean(condRes$pvalue <= 0.05, na.rm = TRUE)

## order genes by significance
condRes <- condRes %>% rownames_to_column('gene') %>% arrange(desc(waldStat))
plot(seq(G), condRes$waldStat)

oo <- order(condRes$waldStat, decreasing = TRUE)
plotSmoothers(sce, assays(sce)$counts,
              gene = rownames(assays(sce)$counts)[oo[450]],
              alpha = 1, border = TRUE)

res.tradeseq <- condRes %>% arrange(desc(waldStat)) %>% 
  mutate(Rank = row_number(), 
         neglogp = -log(pvalue + 1e-16),
         fdr = p.adjust(pvalue, method = 'BH')) %>% as_tibble()

(p4 <- get_ranked_dge_plot(res.tradeseq, method = 'tradeSeq', markers = effect_genes))
res.tradeseq %>% filter(fdr < .05) %>% with(table(gene %in% null_genes)) 

write_rds(res.tradeseq, file = 'output/simulations/tradeSeq_res.rds')
#===================================================================
## run differential testing with Lamian
#===================================================================
## get normalized counts as input 
expr <- assay(vsd)
range(apply(expr, 1, sd))
## pseudotime values
all(bin_meta$sample_bin == colnames(vsd))
pt <- bin_meta$coord_x; names(pt) <- bin_meta$sample_bin
sample_meta <- bin_meta %>% select(condition, sample) %>% distinct()
design <- model.matrix(~ condition, data = sample_meta)
rownames(design) <- sample_meta$sample
testvar <- 2

cell_anno <- data.frame(
  Cell = bin_meta$sample_bin,
  Sample = as.character(bin_meta$sample)
)

# get results for contrast controlling for donor effects
res <- lamian_test(expr = as.matrix(assay(vsd)),
                   cellanno = cell_anno,
                   pseudotime = pt,
                   design = design, 
                   test.type= 'variable', 
                   testvar = testvar, 
                   test.method = 'chisq',
                   permuiter = 100, 
                   overall.only = T,
                   ncores = 1) 


## get differential dynamic genes statistics
stat <- res$statistics %>% rownames_to_column('gene') %>% 
  as_tibble() %>% 
  mutate(neglogp = -log(pval.chisq.overall)) %>% 
  arrange(desc(neglogp)) 


(p5 <- get_ranked_dge_plot(stat, method = 'Lamian', markers = effect_genes))

stat %>% filter(fdr.chisq.overall < .05)
stat %>% filter(fdr.chisq.overall < .05) %>% with(table(gene %in% null_genes)) 

## save results
write_rds(res, 'output/simulations/lamian_DGE_res.rds')

wrap_plots(p4 | p5)

write_rds(dds, 'output/simulations/dds.rds')
