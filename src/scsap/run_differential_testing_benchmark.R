#' read anndata obs and counts, and use that to prepare inputs for analysis at different resolutions
#' tradeSeq is not scalable, I will run it at the resolution of the provided binned isodepths from MultiGaston

library(anndata)
ad <- read_h5ad('~/Desktop/integrated_anndata.h5ad')
obs <-  ad$obs %>% rownames_to_column('spot_id') %>% as_tibble() 
# obs <- read_rds('data/temp/ad_obs.rds')

zeallot::`%<-%`(c(cell_level_meta, bin_level_meta, sample_level_meta),  prepare_metadata(obs))
                  
# cluster_meta, tradeseq_meta),
# (clusters <- unique(cluster_meta$tissue))
# nrow(tradeseq_meta)

# I'll make 24wk the base level for sample meta
# sample_level_meta <- sample_level_meta %>% mutate(age = fct_relevel(age, '24wk'))

## I'll run argo also directly at the scale of the binned isodepth values

## prepare the counts data for running each method
cell_level_counts <- ad$layers['counts']
bin_level_counts <- rowsum(as.matrix(cell_level_counts), group = cell_level_meta$slide_bin) # don't really need this right now
# cluster_level_counts <- rowsum(as.matrix(cell_level_counts), group = cell_level_meta$sample_cluster)

## prepare counts for tradeseq
# tradeseq_counts <- rowsum(as.matrix(cell_level_counts), group = cell_level_meta$slide_tradeseq) 

## filter bin counts to remove sparse bins (that have fewer than 25 cells)
bin_level_cdounts <- bin_level_counts[bin_level_meta$slide_bin, ]

## assert that sample IDs match
stopifnot(all(rownames(bin_level_counts) == bin_level_meta$slide_bin))
stopifnot(all(rownames(cluster_level_counts) == cluster_meta$sample_cluster))
stopifnot(all(rownames(tradeseq_counts) == tradeseq_meta$slide_tradeseq))


## prepare input for running argo
bin_meta <- bin_level_meta %>%
  bind_cols(scale(t(varianceStabilizingTransformation(t(bin_level_counts))), center = T, scale = F))
boxplot(bin_meta$C4b~ bin_meta$condition_bin)

## run with existing pipeline from here on out

cluster_meta <- cluster_meta %>% mutate(condition_cluster = make.names(condition_cluster))
table(cluster_meta$condition_cluster)

diff.testing.methods <- list(
  DESeq2 = run_deseq,
  edgeR = run_edger, 
  Limma = run_limma, 
  tradeSeq = run_tradeseq,
  Lamian = run_lamian
)

## capture the current environment to pass on to differential testing functions
env <- environment() 
## call differential testing passing the local environment
res.list <- map(diff.testing.methods, ~ .x(env))


## append argo results for the same contrst 
KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'
res.argo <- vroom(glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv'))
res.argo <- res.argo %>%
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())

res.list[['ARGO']] <- res.argo

sapply(res.list, function(df) sum(df$fdr < 0.05))
write_rds(res.list, 'output/DGE_testing_res_list.rds')
res.list <- read_rds('output/DGE_testing_res_list.rds')
#===============================================================
## benchmarking
#===============================================================
require(readxl)
require(UpSetR)
require(VennDiagram)
library(PRROC)

## summarize results with upset-plot
res.sig <- lapply(res.list, function(res) filter(res, fdr < 0.01) %>% pull(gene))
sapply(res.sig, length)

pdf('figs/upsetPlot_FDR_01.pdf', width = 8, height = 6)
print(upset(fromList(res.sig), order.by = "freq", nsets = length(res.sig), nintersects = 40), newpage= F)
dev.off()

#===============================================================
## prepare genesets
#===============================================================
# bulk rnaseq dataset
gt.hahn <- read_xlsx('data/sheets/NIHMS1924184-supplement-11.xlsx', sheet = 2, skip = 1) 
gt.hahn <- na.omit(unique(gt.hahn$`Gene Symbol`))

## early single cell study
gt.ximerakis <- read_xlsx('data/sheets/41593_2019_491_MOESM10_ESM.xlsx', sheet = 1)
gt.ximerakis <- na.omit(unique(gt.ximerakis$Gene))

## recent single cell study
gt.jin <- read_xlsx('data/sheets/41586_2024_8350_MOESM4_ESM.xlsx', sheet = 3) 
gt.jin <- na.omit(unique(gt.jin$gene))

gt.sziraki <- read_xlsx('data/sheets/Sziraki.xlsx', sheet = 9, skip = 1)
gt.sziraki <- na.omit(unique(gt.sziraki$`Gene name`))

dataset <- list(
                'Hahn et al' = gt.hahn,
                'Jin et al' = gt.jin,
                'Ximerakis et al' = gt.ximerakis,
                'Sziraki et al' = gt.sziraki
)

sapply(dataset, length)

venn_object <- venn.diagram(
  x = dataset,
  filename = NULL,
  # Styling the circles
  fill = c("#8dd3c7", "#ffffb3", "#bebada", "dodgerblue3"),
  alpha = rep(0.5, 4), # Transparency
  lty = "blank",            # Removes the black borders
  # Styling text
  cex = 1.5,                # Numbers size
  cat.cex = 1.5,             # Category labels size
  margin = 0.15,  # Generates extra whitespace around the diagram (default is usually 0.05)
)

# 4. Plot the image
pdf('figs/benchmark-genesets.pdf', width = 6, height=6)
grid.draw(venn_object)
dev.off()

#============================================================
## get AUPRC curves
#============================================================
## find what genes in our filtered list are not in any of the benchmarks
# markers <- Reduce(union, dataset)
# markers <- intersect(markers, res.list$DESeq2$gene)
# genes <- setdiff(res.sig$ARGO, markers)
# res.argo %>% arrange(desc(neglogp)) %>% filter(gene %in% genes & abs(deviation) > 1)


dataset <- lapply(dataset, function(x) intersect(x, funGP.df$gene))
sapply(dataset, length)

STATISTIC <- 'PR'

compute_ranked_auprc <- function(ranked_list, signature, STATISTIC){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  pg_score <- scores[binary_labels]
  ng_score <- scores[!binary_labels]
  if(STATISTIC == 'PR'){
    pr_obj <- pr.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)  
  }else if(STATISTIC == 'ROC'){
    pr_obj <- roc.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)
  }
  return(pr_obj)
}


## loop through genesets and gather AUPRC curves across methods
for(study in names(dataset)){
  message(study)
  markers <- dataset[[study]]
  pr.curves <- lapply(res.list, function(df)  {
    compute_ranked_auprc(df %>% arrange(desc(neglogp)) %>% pull(gene),
                         signature = markers, 
                         STATISTIC = STATISTIC)
    }
  )
  
  cols <- c( "#bebada", "#8dd3c7", 'lightskyblue2', 'lightgreen', 'lightpink2', 'lightsalmon', 'orchid2')
  
  pdf(glue('figs/{STATISTIC}_{study}.pdf'), width = 6, height = 6)
  par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
  plot(pr.curves[['ARGO']], col = adjustcolor(cols[1], alpha.f = 1), 
       main = glue("{study} \n Precision-Recall Curve"), 
       auc.main = F, lty=2)
  ## plot other curves
  for(ii in seq(1,2)){
    method = names(pr.curves)[ii]
    plot(pr.curves[[method]], col = adjustcolor(cols[ii + 1], alpha.f = .5), lty = 2, lwd = 2, add = TRUE)
  }
  # make legends
  if(STATISTIC == 'PR'){
    lgx <- paste("ARGO (AUPRC = ", round(pr.curves[['ARGO']]$auc.integral, 3), ")")
  }else if(STATISTIC == 'ROC'){
    lgx <- paste("ARGO (AUROC = ", round(pr.curves[['ARGO']]$auc, 3), ")")
  }
  
  for(method in names(pr.curves)[1:5]){
    if(STATISTIC == 'PR'){
      lgx <- c(lgx, 
               glue("{method} (AUPRC = {round(pr.curves[[method]]$auc.integral, 3)}) ")
      )
    }else if(STATISTIC == 'ROC'){
      lgx <- c(lgx, 
               glue("{method} (AUROC = {round(pr.curves[[method]]$auc, 3)}) ")
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
  
}




