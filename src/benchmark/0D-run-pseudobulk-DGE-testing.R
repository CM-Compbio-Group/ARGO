## run DGE testing on pseudo-bulk counts data with different tools and compare performance
# rm(list=ls())
source('src/argo/00-load-packages.R')
CONTRAST <- 'old_vs_adult'

## scripts to gather DGE testing results
source('src/R/get_cct_neglogp.R')
source('src/R/summarize_dge_results.R')
source('src/R/get_ranked_dge_plot.R')

## load DESeq object
library(DESeq2)
domain_dds <- read_rds('data/temp/domain_dds.rds')
# levels(domain_dds$tissue) <- str_replace_all(levels(domain_dds$tissue),  c(" " = '_',  "/" = '_'))
# domain_dds$age_domain <-  factor(paste(domain_dds$age, domain_dds$tissue, sep = '_'))
# write_rds(domain_dds, 'data/temp/domain_dds.rds')
#------------------------------------------------------
## run tests with deSeq
#------------------------------------------------------
design(domain_dds) <- ~ 0 + age_domain

domain_dds <- DESeq(domain_dds)
res <- results(domain_dds)
summary(res, alpha = .05)

## gather results across all bins
res.list <- list()
for(domain in levels(domain_dds$tissue)){
  res.list[[domain]] <- results(domain_dds, contrast= c('age_domain', paste('90wk', domain, sep='_'),
                                                        paste('24wk', domain, sep='_')))
}


head(res.list[[1]])
deseq_res <- summarize_dge_results(res.list, genes = rownames(domain_dds), pval = 'pvalue', lfc = 'log2FoldChange')
write_rds(deseq_res, 'output/deseq_DGE_res.rds')

#----------------------------------
#----------------------------------
## summarize findings
p1 <- get_ranked_dge_plot(deseq_res$dge, method = 'DESeq2')
with(deseq_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
deseq_res$dge %>% filter(neglogp > -log(.05)) %>% with(table(max.lfc > 0)) ## 102 upreg; 94 downreg
#------------------------------------------------------
## run tests with edgeR
#------------------------------------------------------
library(edgeR)

meta <- as_tibble(colData(domain_dds))
y <- DGEList(counts = counts(domain_dds, normalized = F), 
             samples = meta)

## normalize the counts
y <- calcNormFactors(y)
## design matrix for contrasts
design <- model.matrix(~ 0 + age_domain, data = y$samples)
colnames(design) <- gsub("age_domain", "age_", colnames(design))
print(colnames(design))

## run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

## extract results with contrasts
res.list <- list()
for(domain in unique(y$samples$tissue)){
  message(domain)
  contrast_formula <- glue("age_90wk_{domain} - age_24wk_{domain}")
  con <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  ## run the test and extract results
  qlf <- glmQLFTest(fit, contrast = con)
  ## save results to list
  res.list[[domain]] <- topTags(qlf, n = Inf, sort.by = "none")$table
}

names(res.list[[1]])
edger_res <- summarize_dge_results(res.list, genes = rownames(domain_dds), pval = 'PValue', lfc = 'logFC')
write_rds(edger_res, file = 'output/edgeR_DGE_res.rds')

#---------------------------------------------
#---------------------------------------------
## summarize the findings
p2 <- get_ranked_dge_plot(edger_res$dge, method = 'edgeR')
with(edger_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
edger_res$dge %>% filter(neglogp > -log(.05)) %>% with(table(max.lfc > 0)) ## 115 upreg; 85 downreg

#------------------------------------------------------
## run tests with limmma-voom
#------------------------------------------------------
library(limma)

## apply voom-transformation to generate the precision weights
v <- voom(y, design = design, plot = T)
# fit initial linear model
fit <- lmFit(v, design)

## extract results with contrasts
res.list <- list()
for(domain in unique(y$samples$tissue)){
  message(domain)
  contrast_formula <- glue("age_90wk_{domain} - age_24wk_{domain}")
  con <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  ## fit contrast to the linear model
  fit_contrast <- contrasts.fit(fit, contrasts =con)

  ## run empirical bayes moderation
  fit_ebayes <- eBayes(fit_contrast)
  
  ## save results to list
  res.list[[domain]] <- topTable(fit_ebayes, coef = 1, n = Inf, sort.by = "none")
}


print(names(res.list[[1]]))
limma_res <- summarize_dge_results(res.list, rownames(domain_dds), pval = 'P.Value', lfc = 'logFC')
write_rds(limma_res, file = 'output/limma_DGE_res.rds')
#---------------------------------------------
#---------------------------------------------
## summarize the findings
p3 <- get_ranked_dge_plot(limma_res$dge, method = 'limma-voom')
with(limma_res$dge, plot(neglogp~max.lfc, pch=20, col = adjustcolor("red", alpha.f = 0.5)))
limma_res$dge %>% filter(neglogp > -log(.05)) %>% with(table(max.lfc > 0)) ## 103 upreg; 87 downreg
#---------------------------------------------
#---------------------------------------------

pdf('figs/diff_expression_testing.pdf', width = 12)
plot(wrap_plots(p1 | p2 | p3))
dev.off()
