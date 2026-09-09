## normalize the curves 
y_shape <- y / apply(abs(y), 1, max)
dist_mat <- dist(y_shape) / sqrt(diff(range(x)))

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


pdf('figs/hclust_shape_silhouette.pdf', width = 10, height = 4)
par(mfrow = c(1,3), cex.axis = 1, 
    cex.lab = 1, 
    cex.main = 1,
    mar = c(4, 4.1, 2.0, 2.1) )
generate_elbow(hc_comp, 'complete')
generate_elbow(hc_avg, 'average')
generate_elbow(hc_sing, 'single')
dev.off()


table(cutree(hc_comp, k = 4), sig_genes$net_auec > 0)
table(cutree(hc_comp, k = 7), sig_genes$net_auec > 0)

#--------------------------------
#--------------------------------
## fit with pam clustering
fit_pam_y <- pam(y, k = 7, metric = "euclidean")
fit_pam_ys <- pam(y_shape, k = 7, metric = "euclidean")


# 2. Extract cluster assignments and medoid curves
cluster_ids   <- fit_pam_ys$clustering
medoid_curves <- fit_pam_ys$medoids     # The 7 actual representative gene curves
medoid_row_idx <- fit_pam_ys$id.med     # Row indices of the 7 medoid genes

cxid <- fit_pam_y$clustering
table(cxid, cluster_ids)
cid   <- cutree(hc_comp, k = 7)
table(cluster_ids, cid)

# 3. View cluster sizes
table(cluster_ids, sig_genes$net_auec > 0)

## visualize clusters
pdf('figs/k7_shape_PAM_partition.pdf', width = 6, height = 6)
cols <- brewer.pal(7, "Set2")
matplot(x, t(y_shape), type = "l", col = adjustcolor(cols[cluster_ids], alpha.f = 0.3),
        lty = 1, xlab = "Isodepth Axis", ylab = "Expression",
        main = "7-Cluster Medoid Partition on Shape")

# Overlay the 7 actual medoids
matlines(x, t(medoid_curves), type = "l", col = cols, lwd = 5, lty = 2)
legend("topright", 
       legend = paste("Cluster", 1:7), 
       col = cols, 
       lwd = 3, 
       lty = 2, 
       bty = "n",         # Removes the box border around the legend
       cex = 0.8)         # Controls the text size
dev.off()


sig_genes$shape_cluster <- cluster_ids

pdf('figs/shape_cluster_summary.pdf', width = 8, height = 4)
par(mfrow = c(1,2), cex.axis = 1, cex.lab = 1, cex.main = 1,
    mar = c(4, 4.1, 2.0, 2.1) )
plot(sig_genes$net_auec~ sig_genes$shape_cluster, xlab = 'Cluster', ylab = 'Net AUEC', pch = 18, cex=1.5)
plot(sig_genes$index~ sig_genes$shape_cluster, xlab = 'Cluster',
     ylab = 'Index of max deviation', pch = 18, cex=1.5)
dev.off()