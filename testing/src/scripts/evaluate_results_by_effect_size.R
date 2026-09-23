#===================================================================
## compute eval stats by strength of signals
#===================================================================

## gather the metrics



eval.list <- lapply(dge.res.list, compute_eval_metrics, dplyr::filter(effect_info, amplitude == 1) %>% pull(gene), null_genes)

## gather the metrics
print(eval.metrics <- do.call('rbind', lapply(eval.list, '[[', 'metrics')))
pr.curves <- lapply(eval.list, '[[', 'probj')

SX <- glue('pancreatic_endocrinogenesis-{iter}')
cols <- c( 'burlywood2', 'lightskyblue2', 'lightgreen', 'lightpink', 'lightsalmon', 'plum3', 'orchid3')

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

