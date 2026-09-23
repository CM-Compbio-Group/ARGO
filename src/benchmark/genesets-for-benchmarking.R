## find no. of genes in each benchmarking dataset
source('src/scripts/00-load-packages.R')
require(VennDiagram)
library(PRROC)

CONTRAST <- 'old_vs_adult'
KERNEL <- 'RBF'
res.list <- read_rds(glue('output/{CONTRAST}_DGE_res_list.rds'))
res.fungp <- res.list$funGP
res.lamian <- res.list$lamian

#===================================================================
#===================================================================
## summarize results with upset-plot
require(UpSetR)
res.sig <- lapply(res.list, function(res) filter(res, neglogp > -log(0.01)) %>% pull(gene))
sapply(res.sig, length)

pdf('figs/upsetPlot_pval_01.pdf', width = 8, height = 6)
print(upset(fromList(res.sig), order.by = "freq"), newpage= F)
dev.off()

#===================================================================
#===================================================================
## compare rankings of the top funGP-ranked genes across other methods
genes <- res.fungp$gene[1:10]
lapply(res.list, function(res) filter(res, gene %in% genes))


#===================================================================
#===================================================================

# gt.hahn.age_corr <- read_xlsx('data/sheets/NIHMS1924184-supplement-11.xlsx', sheet = 3, skip = 1) 
# gt.hahn.age_corr <- gt.hahn.age_corr %>% filter(`Gene Symbol` %in% res.fungp$gene)
# gt.hahn.age_corr <- na.omit(unique(gt.hahn.age_corr$`Gene Symbol`))

# bulk rnaseq dataset
gt.hahn.age_diff <- read_xlsx('data/sheets/NIHMS1924184-supplement-11.xlsx', sheet = 2, skip = 1) 
# gt.hahn.age_diff <- gt.hahn.age_diff %>% filter(`Gene Symbol` %in% res.fungp$gene)
gt.hahn.age_diff <- na.omit(unique(gt.hahn.age_diff$`Gene Symbol`))

## early single cell study
gt.ximerakis <- read_xlsx('data/sheets/41593_2019_491_MOESM10_ESM.xlsx', sheet = 1)
# gt.ximerakis <- gt.ximerakis %>% filter(Gene %in% res.fungp$gene)
gt.ximerakis <- na.omit(unique(gt.ximerakis$Gene))

## recent single cell study
gt.jin <- read_xlsx('data/sheets/41586_2024_8350_MOESM4_ESM.xlsx', sheet = 3) 
# gt.jin <- gt.jin %>% filter(gene %in% res.fungp$gene)
gt.jin <- na.omit(unique(gt.jin$gene))


## visium ST study
gt.wang.1 <- read_xlsx('data/sheets/wang_et_al_DARs.xlsx', sheet = 1)
gt.wang.2 <- read_xlsx('data/sheets/wang_et_al_DARs.xlsx', sheet = 2)
gt.wang <- na.omit(unique(c(gt.wang.1$`gene ID`, gt.wang.2$`gene ID`)))
# gt.wang <- intersect(gt.wang, res.fungp$gene)


## IRISeq study
gt.abdulraouf <- read_xlsx('data/sheets/IRISeq.xlsx', sheet = 5, skip=1)
# gt.abdulraouf <- filter(gt.abdulraouf, gene_id %in% res.fungp$gene)
gt.abdulraouf <- na.omit(unique(gt.abdulraouf$gene_id))

dataset <- list('hahn' = gt.hahn.age_diff, 
                'jin' = gt.jin,
                'ximerakis' = gt.ximerakis,
                'wang' = gt.wang,
                'abdulraouf' = gt.abdulraouf
                )

sapply(dataset, length)

#===================================================================
#===================================================================

venn_object <- venn.diagram(
  x = dataset,
  category.names = c("Hahn et al", "Jin et al", "Ximerakis et al", 
                     'Wang et al', 'Abdulraouf et al'),
  filename = NULL,
  
  # Styling the circles
  fill = c("#8dd3c7", "#ffffb3", "#bebada", 'lightpink', 'orchid2'),
  alpha = rep(0.5, 5), # Transparency
  lty = "blank",            # Removes the black borders
  
  # Styling text
  cex = 1.5,                # Numbers size
  cat.cex = 1.5,             # Category labels size
  
  margin = 0.15,  # Generates extra whitespace around the diagram (default is usually 0.05)
  cat.pos = c(0, 280, 240, 120, 40), # Angles (in degrees) to position the 5 labels around the circles
  cat.dist = c(0.25, 0.25, 0.25, 0.25, 0.25) # Distance of labels from the circle edges
)

# 4. Plot the image
pdf('figs/benchmark-genesets-complete.pdf', width = 8, height=8)
grid.draw(venn_object)
dev.off()

#===================================================================
#===================================================================
source('src/R/get_ranked_dge_plot.R')

markers <- Reduce(intersect, dataset)
markers <- intersect(markers, res.fungp$gene)

## find what genes in our filtered list are not in any of the benchmarks
markers <- Reduce(union, dataset)
markers <- intersect(markers, res.fungp$gene)
genes <- setdiff(res.sig$funGP, markers)
res.fungp %>% arrange(desc(neglogp)) %>% filter(gene %in% genes & abs(deviation) > 1)


#===================================================================
#===================================================================
SX <- 'intersection'
statistic <- 'PR' ## whether to plot AUPRC or AUROC curves. Set to PR or ROC
markers <- dataset[['abdulraouf']]

p1 <- get_ranked_dge_plot(res.list$funGP, method = 'funGP', markers = markers)
p2 <- get_ranked_dge_plot(res.list$lamian, method = 'lamian', markers = markers)
pdf(glue('figs/{SX}_ranked-statistics.pdf'), width = 12)
wrap_plots(p1 | p2)
dev.off()

#===================================================================
compute_ranked_auprc <- function(ranked_list, signature, statistic){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  pg_score <- scores[binary_labels]
  ng_score <- scores[!binary_labels]
  if(statistic == 'PR'){
    pr_obj <- pr.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)  
  }else if(statistic == 'ROC'){
    pr_obj <- roc.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)
  }
  return(pr_obj)
}


pr.curves <- lapply(res.list, function(df)  {
  compute_ranked_auprc(df %>% arrange(desc(neglogp)) %>% pull(gene), signature = markers, statistic = statistic)
  }
)


cols <- c( 'red3', 'lightskyblue2', 'lightgreen', 'lightpink2', 'lightsalmon', 'orchid2')

pdf(glue('figs/{CONTRAST}_{statistic}_comparative_{SX}.pdf'), width = 8, height = 8)
par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
plot(pr.curves[['funGP']], col = adjustcolor(cols[1], alpha.f = 1), main = glue("{SX} \n Precision-Recall Curve"), 
     auc.main = F, lty=2)
## plot other curves
for(ii in seq(1,5)){
  method = names(pr.curves)[ii]
  plot(pr.curves[[method]], col = adjustcolor(cols[ii + 1], alpha.f = .5), lty = 2, lwd = 2, add = TRUE)
}
## make legends
if(statistic == 'PR'){
  lgx <- paste("FunGP (AUPRC = ", round(pr.curves[['funGP']]$auc.integral, 3), ")")  
}else if(statistic == 'ROC'){
  lgx <- paste("FunGP (AUPRC = ", round(pr.curves[['funGP']]$auc, 3), ")")
}

for(method in names(pr.curves)[1:5]){
  if(statistic == 'PR'){
    lgx <- c(lgx, 
             glue("{method} (AUPRC = {round(pr.curves[[method]]$auc.integral, 3)}) ")
    )
  }else if(statistic == 'ROC'){
    lgx <- c(lgx, 
             glue("{method} (AUPRC = {round(pr.curves[[method]]$auc, 3)}) ")
    )
  }
  
  
}
legend("topright", cex=.8,
       legend = lgx,
       col = cols, 
       lty = 2,
       lwd = 2,             # Line width
       bty = "n")           # Removes the box border around the legend (optional)

dev.off()
