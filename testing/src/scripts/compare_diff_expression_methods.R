## gather all results and make plots
library(readxl)
library(PRROC)
library(extrafont)

## fungp experiment variants
EXPT1 <- 'single_noise'
EXPT2 <- 'conditional_noise'

# res.deseq <- read_rds('output/simulations/deseq_DGE_res.rds')[['dge']]
# res.edger <- read_rds('output/simulations/edgeR_DGE_res.rds')[['dge']]
# res.limma <- read_rds('output/simulations/limma_DGE_res.rds')[['dge']]
# res.tradeseq <- read_rds('output/simulations/tradeSeq_res.rds')
# res.lamian <- read_rds('output/simulations/lamian_DGE_res.rds')
# res.fungp.expt1 <- read_rds(glue('output/simulations/funGP_res_{EXPT1}.rds')) 
# res.fungp.expt2 <- read_rds(glue('output/simulations/funGP_res_{EXPT2}.rds')) 


## extract statistics for lamian
res.lamian <- res.lamian$statistics %>% 
  rownames_to_column('gene') %>% 
  as_tibble() %>% 
  mutate(neglogp = -log(pval.chisq.overall)) %>% 
  rename(fdr = fdr.chisq.overall) %>% 
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number())


res.list <- list('deseq' = res.deseq,
                 'edgeR' = res.edger,
                 'limma' = res.limma,
                 'tradeseq' = res.tradeseq,
                 'lamian' = res.lamian,
                 'funGP_EXPT1' = res.fungp.expt1,
                 'funGP_EXPT2' = res.fungp.expt2)

fdr_thresh <- 0.05
sapply(res.list, function(df) sum(df[['fdr']] < fdr_thresh))

#====================================================#
## compute the auprc curves
#====================================================#
compute_ranked_auprc <- function(ranked_list, signature){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  positive_score <- scores[binary_labels]
  negative_score <- scores[!binary_labels]
  pr_obj <- pr.curve(scores.class0 = positive_score, scores.class1 = negative_score, curve = T)
  return(pr_obj)
}


markers <- effect_genes
SX <- 'scDesign3 pancreatic endocrinogenesis'
cols <- c( 'burlywood2', 'lightskyblue2', 'lightgreen', 'lightpink', 'lightsalmon', 'plum3', 'orchid3')

## I only want to rank by p-values
pr.curves <- lapply(res.list, function(df)  {
    compute_ranked_auprc(df %>% arrange(desc(neglogp)) %>% pull(gene), signature = markers)
  }
)

names(pr.curves)[c(6,7)] <- c('funGP_single_noise', 'funGP_conditional_noise')

pdf(glue('figs/AUPRC_{SX}.pdf'), width = 8, height = 8)
par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
plot(pr.curves[['deseq']], col = adjustcolor(cols[1], alpha.f = 1), 
     main = glue("{SX} \n Precision-Recall Curve"), 
     auc.main = F, lty=2, lwd = 2)
## plot other curves
for(ii in seq(2,7)){
  method = names(pr.curves)[ii]
  plot(pr.curves[[method]], col = adjustcolor(cols[ii + 1], alpha.f = .5), lty = 2, lwd = 2, add = TRUE)
}

## make legends
lgx <- glue("DESeq2 (AUPRC = ", round(pr.curves[['deseq']]$auc.integral, 3), ")")
for(method in names(pr.curves)[2:7]){
  lgx <- c(lgx, 
              glue("{method} (AUPRC = {round(pr.curves[[method]]$auc.integral, 3)}) ")
              )
}
legend("bottomleft", cex=.8,
       legend = lgx,
       col = cols, 
       lty = 2,
       lwd = 1,             # Line width
       bty = "n")           # Removes the box border around the legend (optional)

dev.off()


#====================================================#
## compute table of recovery rates for each method
#====================================================#
## compute performance metrics for each method
compute_recovery_metrics <- function(df, markers){
  
  called <- filter(df, fdr < 0.05) %>% pull(gene)
  not_called <- filter(df, fdr >= 0.05) %>% pull(gene)
  TP = length(intersect(called, markers))
  FP = length(setdiff(called, markers))
  TN = length(setdiff(not_called, markers))
  FN = length(intersect(not_called, markers))
  
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  false_positive_rate = FP/ (FP + TN)
  false_negative_rate = FN/ (FN + TP)
  
  precision = TP / (TP + FP)
  empirical_FDR = FP / (TP + FP)
  
  accuracy = (TP + TN) /  (TP + FP + TN + FN)
  F1 = (2 * TP) /  (2 * TP + FP + FN)
  
  
  return(c('TP' = TP, 'FP' = FP, 'TN' = TN, 'FN' = FN, 
           'sensitivity' = sensitivity, 
           'specificity' = specificity, 
           'false_positive_rate' = false_positive_rate,
           'false_negative_rate' = false_negative_rate,
           'precision' = precision,
           'empirical_FDR' = empirical_FDR,
           'accuracy' = accuracy,
           'F1' = F1))
}

metrics <- lapply(res.list, compute_recovery_metrics, markers)
metrics <- do.call('rbind', metrics)
write.csv(metrics, file = glue('output/simulations/{SX}_performance_metrics.csv'), quote = F)
