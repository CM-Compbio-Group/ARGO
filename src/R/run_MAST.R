#---------------------------------------------
## run the tests using MAST
#---------------------------------------------
library(MAST)
## extract log-cpm counts
cpm_counts <- cpm(counts(domain_dds, normalized = F), log = F)
log_cpm <- log2(cpm_counts + 1)

## construct single cell assay
sca <- FromMatrix(exprsArray = log_cpm, cData = data.frame(colData(domain_dds)))

# set the cellular detection rate for each gene
colData(sca)$cdr <- colMeans(assay(sca) > 0) # proportion of detected genes per sample

## fit the MAST hurdle model
fit <- zlm(~ 0 + age_domain + cdr, sca)

## clean up variable names in fit
colnames(fit@coefC) <- gsub('age_domain', 'age_', colnames(fit@coefC))
colnames(fit@coefD) <- gsub('age_domain', 'age_', colnames(fit@coefD))

## extract results looping over contrasts
res.list <- list()
for(domain in unique(colData(sca)$tissue)){
  message(domain)
  baseline_coef  <- glue("age_24wk_{domain}")
  treatment_coef <- glue("age_90wk_{domain}")
  contrast_formula <- glue("age_90wk_{domain} - age_24wk_{domain}")
  
  ## evaluate LRT on the contrast
  lrt <- lrTest(fit, Hypothesis(contrast_formula))
  
  # #. Extract the combined Hurdle model statistics
  # 'hurdle' contains the combined p-values from both components
  lrt_df <- as.data.frame(lrt[, "hurdle", ])
  lrt_df$primerid <- rownames(lrt_df)
  
  # #. Get Log-Fold Changes (LogFC)
  # MAST requires getLogFC to compute explicit contrast-level fold changes
  logfc_res <- getLogFC(fit, contrast0 = Hypothesis(baseline_coef), contrast1 = Hypothesis(treatment_coef))
  
  # #. Merge the Test Statistics and LogFC values into a cohesive table
  merged_res <- merge(logfc_res, lrt_df, by = "primerid")
  
  # Clean up column names and compute an Adjusted P-Value (FDR)
  colnames(merged_res)[colnames(merged_res) == "Pr(>Chisq)"] <- "p.value"
  merged_res$fdr <- p.adjust(merged_res$p.value, method = "BH")
  
  # #. Save results to the list
  res.list[[domain]] <- merged_res
  
}
