## cluster the significant results and perform downstream inference
# source('src/scripts/00-load-packages.R')
# library(fdacluster)
# library(future)
library(patchwork)

# pvalsDf <- vroom('data/processed/funGP_testing_old_vs_juvenile.csv')
# sig_genes <- vroom('data/processed/sig_de_genes.csv')

#--------------------------------
#--------------------------------
## filter genes to fit
sig_genes <- filter(pvalsDf,(neglogp > -log(.05)) & abs(net_auec) > .25)
table(sig_genes$net_auec > 0)

## get the curves for these genes
y_curves <- do.call('rbind',lapply(sig_genes$gene, function(gene) exact_pvals[[gene]][['statistic']]))
x <- x_test

#--------------------------------
## find optimal no. of clusters w/out warping using hclust
#--------------------------------
# 1. Compute normalized L2 distance matrix once
# (L2 distance divided by sqrt(range))
dist_mat <- dist(y_curves) / sqrt(diff(range(x)))
dim(dist_mat)

# 2. Fit the hierarchical models
hc_comp <- hclust(dist_mat, method = "complete")
hc_avg  <- hclust(dist_mat, method = "average")
hc_sing <- hclust(dist_mat, method = "single")

# 3. Quickly compute average silhouettes across k = 2:15
library(cluster)
k_range <- 2:15

generate_elbow <- function(hc, method){
  sil_results <- sapply(k_range, function(k) {
    message(k)
    sil <- silhouette(cutree(hc, k = k), dist_mat)
    mean(sil[,'sil_width'])
  })
  
  plot(k_range, sil_results, main = glue('{method} clustering'),
       pch=20,cex=1.5, xlab = 'no. of clusters', ylab = 'mean silhouette width')
}


pdf('figs/hclust_silhouette.pdf', width = 10, height = 4)
par(mfrow = c(1,3), cex.axis = 1, 
    cex.lab = 1, 
    cex.main = 1,
    mar = c(4, 4.1, 2.0, 2.1) )
generate_elbow(hc_comp, 'complete')
generate_elbow(hc_avg, 'average')
generate_elbow(hc_sing, 'single')
dev.off()

## compare expression trends by cluster
table(cutree(hc_comp, k=7), sig_genes$net_auec > 0)
#--------------------------------
#--------------------------------
## do the same with the fdacluster package
ncores <- max(parallel::detectCores() - 4L, 1L)
plan(multisession, workers = ncores)
start_time <- Sys.time()
amp_mcaps <- compare_caps(
  x = x,
  y_curves = y_curves,
  n_clusters = 2:10,
  metric = "l2", 
  warping_class = 'none',
  # warping_class = c("none","shift"), #, "dilation", "affine"),
  # clustering_method = 'hclust-complete',
  clustering_method = 'kmeans',
  centroid_type = "medoid",
  cluster_on_phase = F
)
end_time <- Sys.time()
end_time - start_time
plan(sequential)

## silhouettes not computed by default with kmeans clustering
# 2. Extract silhouette scores for each k model
sil_list <- lapply(amp_mcaps$caps_obj, function(m) {
  sil <- silhouette(m$memberships, dist_mat)
  data.frame(
    k = as.factor(length(unique(m$memberships))),
    sil_width = sil[, "sil_width"]
  )
})
sil_df <- do.call(rbind, sil_list)

# 3. Plot silhouette distribution (Boxplot across k)
pdf('figs/kmeans_silhouette_distribution.pdf', width = 6, height = 4)
boxplot(sil_width ~ k, data = sil_df,
        col = "lightblue", pch = 19,
        xlab = "Number of Clusters (k)",
        ylab = "Silhouette Width",
        main = "Silhouette Distribution Across K-Means Models")
grid()
dev.off()

# 4. Plot mean silhouette trendline
mean_sils <- tapply(sil_df$sil_width, sil_df$k, mean)
pdf('figs/kmeans_silhouette_means.pdf', width = 6, height = 4)
plot(names(mean_sils), mean_sils, type = "b", pch = 19, col = "steelblue", lwd = 2,
     xlab = "Number of Clusters (k)", ylab = "Mean Silhouette Width",
     main = "Mean Silhouette vs K") 
grid()
dev.off()

plot(amp_mcaps, validation_criterion = "wss", what = "distribution")
plot(amp_mcaps, validation_criterion = "wss", what = "mean")



## plot the resulting clusters
out <- amp_mcaps$caps_obj[[6]]
pdf('figs/k8_fits.pdf', width = 8, height = 8)
plot(out, type = 'amplitude')
dev.off()

table(out$memberships, sig_genes$net_auec > 0)

#--------------------------------
#--------------------------------
## check output of fitting with time-warps
out <- fdakmeans(
  x,
  y_curves,
  n_clusters = 7,
  centroid_type = "mean",
  warping_class = "shift",
  metric = "l2",
  cluster_on_phase = T, 
  # linkage_criterion = 'average'
)

plot(out)

#--------------------------------
#--------------------------------
## fit final curve and visualize it
# 1. Run exact medoid clustering (takes < 0.05 seconds)
fit_pam <- pam(y_curves, k = 7, metric = "euclidean")

# 2. Extract cluster assignments and medoid curves
cluster_ids   <- fit_pam$clustering
medoid_curves <- fit_pam$medoids     # The 7 actual representative gene curves
medoid_row_idx <- fit_pam$id.med     # Row indices of the 7 medoid genes

# 3. View cluster sizes
table(cluster_ids, sig_genes$net_auec > 0)

## visualize clusters
pdf('figs/k7_PAM_partition.pdf', width = 8, height = 8)
cols <- rainbow(7)
matplot(x, t(y_curves), type = "l", col = adjustcolor(cols[cluster_ids], alpha.f = 0.3),
        lty = 1, xlab = "Isodepth Axis", ylab = "Expression",
        main = "7-Cluster Medoid Partition (No Warping)")

# Overlay the 7 actual medoids
matlines(x, t(medoid_curves), type = "l", col = "black", lwd = 2, lty = 2)
dev.off()

#--------------------------------
## make inferences with fitted clusters
#--------------------------------
sig_genes$cluster <- cluster_ids

pdf('figs/cluster_summary.pdf', width = 8, height = 4)
par(mfrow = c(1,2), cex.axis = 1, cex.lab = 1, cex.main = 1,
    mar = c(4, 4.1, 2.0, 2.1) )
plot(sig_genes$net_auec~ sig_genes$cluster, xlab = 'Cluster', ylab = 'Net AUEC', pch = 18, cex=1.5)
plot(sig_genes$index~ sig_genes$cluster, xlab = 'Cluster',
     ylab = 'Index of max deviation', pch = 18, cex=1.5)
dev.off()

write_csv(sig_genes, 'data/processed/sig_de_genes.csv')
#--------------------------------
## run go enrichment for clustered-genes
#--------------------------------
library(enrichR)
dir.create('figs/enrichr')

View(listEnrichrDbs())

dbs <- c("GO_Biological_Process_2026", "TISSUES_Curated_2025", "MSigDB_Hallmark_2020",
         "CellMarker_2024", "GTEx_Aging_Signatures_2021", "Tabula_Muris",
         "Reactome_Pathways_2024")

bg <- pvalsDf$gene
for(clust in seq(7)){
  message(clust)
  genes <- filter(sig_genes, shape_cluster == clust) %>% pull(gene)
  enr <- enrichr(genes, dbs)
  
  for(db in dbs){
    if(nrow(enr[[db]])){
      message(db)  
      pdf(glue('figs/enrichr/cluster_{clust}_{db}.pdf'), width=10, height = 6)
      plot(plotEnrich(enr[[db]], showTerms = 5, numChar = 50, 
                 y_curves = "Count", orderBy = "Adjusted.P.value", 
                 title = glue('Cluster-{clust} {db}'))) + 
        theme(axis.text.y_curves = element_text(size = 16)) # Change 12 to your desired font size
      
      dev.off()
    }
  }
  # printEnrich(enr, outFile = "txt", prefix = clust)
}

#--------------------------------
## visualize medoids in space
#--------------------------------
library(Seurat)
dxlist <- read_rds('data/processed/slide_list.rds')

filter(sig_genes, shape_cluster == 2)
sig_genes[medoid_row_idx, ]

medoid_genes <- (sig_genes %>% pull(gene))[medoid_row_idx]
## for these genes, identify rois to plot
idx_medoid <- lapply(exact_pvals[medoid_genes], function(kxx) which(kxx$statistic > kxx$band))
idx_medoid[5:7] <- lapply(exact_pvals[medoid_genes[5:7]], function(kxx) which(kxx$statistic < -kxx$band))

## get corresponding isodepth values
idx_medoid <- lapply(idx_medoid, function(ix) x_test[ix])

## how large are the rois
sapply(idx_medoid, length)


## indicate ROIs by size / alpha
plot_roi <- function(ad, roi, gene){
  mask <- (ad$binned_isodepth %in% roi)
  alpha_map <- ifelse(mask, .25, .05)
  p <- FeaturePlot(ad, gene, alpha = alpha_map, 
              pt.size = alpha_map,
              reduction = 'X_spatial_coords', 
              max.cutoff = 'q9') + 
    labs(x = '', y_curves = '') + ggtitle('')
  return(p)
}

medoid_figs <- lapply(medoid_genes, function(gene){
  lapply(dxlist, plot_roi, idx_medoid[[gene]], gene)
  }
)
 
names(medoid_figs) <- medoid_genes
figs <- lapply(medoid_genes, function(gene) {
  row_plot <- wrap_plots(medoid_figs[[gene]], ncol = 3) + 
    # plot_layout(guides = "collect") + 
    plot_annotation(title = gene)
  # Protect the row assembly from being cleared by the top-level layout
  wrap_elements(row_plot)
  }
)


pdf(glue('figs/{contrast}_upreg_medoids.pdf'), width = 12, height = 15)
p <- wrap_plots(figs[1:4], nrow = 4) + 
  plot_annotation(title = glue('{contrast} upregulated medoids'))
print(p)
dev.off()


pdf(glue('figs/{contrast}_dnreg_medoids.pdf'), width = 12, height = 12)
p <- wrap_plots(figs[5:7], nrow = 3) + 
  plot_annotation(title = glue('{contrast} downregulated medoids'))
print(p)
dev.off()

#--------------------------------
## visualize confidence bands for the medoids
#--------------------------------
conf_band_plots <- lapply(medoid_genes, function(gene) 
  plot_CI_bands(gene, x = x_test, 
                band = exact_pvals[[gene]][['band']],
                stat = exact_pvals[[gene]][['statistic']], 
                alpha = 95)
)


ggsave(wrap_plots(conf_band_plots, ncol = 4) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{contrast}_medoids_conf_bands.pdf')),
        width=18,height=8)

#--------------------------------
## compare w/ deseq results
#--------------------------------
dfs2 <- vroom('data/processed/DESeq_DGE_bin_old_vs_juv.csv') %>%
  column_to_rownames('gene')

dfs2 <- dfs2[sig_genes$gene, ]
all(rownames(dfs2) == sig_genes$gene)

plot(dfs2$rank~ sig_genes$Rank)
plot(dfs2$stats~ sig_genes$stats)
table(dfs2$maxlogFC > 0, sig_genes$net_auec > 0)

gx <- sig_genes$gene[which((dfs2$maxlogFC > 0) & (sig_genes$net_auec < 0))]
dfs2[gx, ]
filter(sig_genes, gene %in% gx)
