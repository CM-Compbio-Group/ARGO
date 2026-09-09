## cluster genes based on shape of the statistic-curve
source('src/scripts/00-load-packages.R')
source('src/R/compute_exceedence_statistics.R')
library(cluster) 
library(RColorBrewer)

exact_pvals <- read_rds(glue('output/funGP_testing_{KERNEL}_{CONTRAST}_res.rds'))
pvalsDf <- vroom(glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv'))
## filter to significant genes
sig_genes <- filter(pvalsDf,(neglogp > -log(.05)) & abs(net_auec) > .25)
table(sig_genes$net_auec > 0)

## get the curves for these genes
y_curves <- do.call('rbind',lapply(sig_genes$gene, function(gene) exact_pvals[[gene]][['statistic']]))
grid_x <- x_test; 

#--------------------------------
## find optimal no. of clusters for curves
#--------------------------------
# 1. Compute the distance matrix (equivalent to unwarped L2)
dist_mat <- dist(y_curves, method = "euclidean")

# 2. Initialize variables to store results
k_range <- 2:15
pam_results <- list()
sil_widths <- numeric(length(k_range))

# 3. Iterate through cluster numbers to find the best fit
for (i in seq_along(k_range)) {
  k <- k_range[i]
  
  # Fit Partitioning Around Medoids (PAM)
  # diss = TRUE tells it we are passing a distance matrix
  fit <- pam(dist_mat, k = k, diss = TRUE)
  
  pam_results[[i]] <- fit
  
  # Extract the average silhouette width (closer to 1 is better)
  sil_widths[i] <- fit$silinfo$avg.width
}

pdf(glue('figs/{KERNEL}_{CONTRAST}_PAM_cluster_silhouette.pdf'), width = 8, height = 6)
plot(k_range, sil_widths, type = "b", pch = 19, 
     xlab = "Number of Clusters (K)", 
     ylab = "Average Silhouette Width",
     main = "K-Medoids Optimization")
dev.off()

#--------------------------------
## fit k-clusters and compare with diff expression direction
#--------------------------------
K = 6
fit_k <- pam(y_curves, k = K, metric = "euclidean")

# 2. Extract cluster assignments and medoid curves
cluster_ids   <- fit_k$clustering
medoid_curves <- fit_k$medoids     # The 7 actual representative gene curves
medoid_row_idx <- fit_k$id.med     # Row indices of the 7 medoid genes

# 3. View cluster sizes
table(cluster_ids, sig_genes$net_auec > 0)

## visualize clusters
pdf(glue('figs/{KERNEL}_{CONTRAST}_PAM_{K}partition.pdf'), width = 8, height = 8)
cols <- brewer.pal(K, "Set2")
matplot(grid_x, t(y_curves), type = "l", 
        col = adjustcolor(cols[cluster_ids], alpha.f = 0.3),
        lty = 1, xlab = "Isodepth Axis", ylab = "Expression",
        main = glue("{K}-Cluster Medoid Partition (No Warping)"))

# Overlay the K actual medoids
matlines(grid_x, t(medoid_curves), type = "l", col = cols, lwd = 3, lty = 1)
dev.off()

#--------------------------------
## fit clusters based on shapes
#--------------------------------
y_shape <- y_curves / apply(abs(y_curves), 1, max)
dist_mat <- dist(y_shape, method = 'euclidean') 

k_range <- 2:15
pam_results <- list()
sil_widths <- numeric(length(k_range))
# 3. Iterate through cluster numbers to find the best fit
for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- pam(dist_mat, k = k, diss = TRUE)
  pam_results[[i]] <- fit
  sil_widths[i] <- fit$silinfo$avg.width
}

pdf(glue('figs/{KERNEL}_{CONTRAST}_PAM_shape_silhouette_auc_015.pdf'), width = 8, height = 6)
plot(k_range, sil_widths, type = "b", pch = 19, 
     xlab = "Number of Clusters (K)", 
     ylab = "Average Silhouette Width",
     main = "K-Medoids Optimization \n on shape of statistic")
dev.off()

## fit to K curves
K = 13
fit_k <- pam(y_shape, k = K, metric = 'euclidean')
cluster_ids   <- fit_k$clustering
table(cluster_ids, sig_genes$net_auec > 0)

medoid_curves <- fit_k$medoids     # The K actual representative gene curves
medoid_row_idx <- fit_k$id.med     # Row indices of the K medoid genes
table(cluster_ids, sig_genes$net_auec > 0)
## visualize clusters
pdf(glue('figs/{KERNEL}_{CONTRAST}_PAM_{K}partition_on_shape_auc_025.pdf'), width = 8, height = 8)
cols <- brewer.pal(K, "Set3")
matplot(grid_x, t(y_shape), type = "l", 
        col = adjustcolor(cols[cluster_ids], alpha.f = 0.3),
        lty = 1, xlab = "Isodepth Axis", ylab = "Expression",
        main = glue("{K}-Cluster Medoid Partition (No Warping)"))

# Overlay the K actual medoids
matlines(grid_x, t(medoid_curves), type = "l", col = cols, lwd = 3, lty = 1)
legend("topright", 
       legend = paste("Cluster", 1:K), col = cols, lwd = 3, lty = 2, 
       bty = "n", cex = 0.8)  
              
dev.off()

#--------------------------------
## visualize clusters in heatmap
#--------------------------------
cluster_assignments <- factor(fit_k$clustering)
names(cluster_assignments) <- sig_genes$gene
cols <- setNames(brewer.pal(K, "Set3"), 1:K)

row_anno <- rowAnnotation(
  PAM = cluster_assignments,
  col = list(PAM = cols), 
  annotation_name_side = "top"
)

pdf(glue('figs/{KERNEL}_{CONTRAST}_binary_heatmap_clustered.pdf'), width = 12, height = 12)
Heatmap(stat_binary[sig_genes$gene, ], right_annotation = row_anno,
        cluster_columns = F, cluster_rows = T,
        row_split = cluster_assignments,
        cluster_row_slices = F,
        name = 'binary \nstatistic',
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=4))
dev.off()

myc <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue3", "lightgray", "red3"))
pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_heatmap_clustered_{K}.pdf'), width = 10, height = 10)
Heatmap(stat_lfc[sig_genes$gene, ], right_annotation = row_anno,
        cluster_columns = F, cluster_rows = T,
        row_split = cluster_assignments,
        col = myc,
        cluster_row_slices = F,
        name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=4))
dev.off()

sig_genes$cluster <- cluster_assignments
write_csv(sig_genes, file = glue('data/processed/{KERNEL}_{CONTRAST}_sigGenes.csv'))

#--------------------------------
## run enrichment analysis with gene clusters
#--------------------------------
library(enrichR)
# dir.create('figs/enrichr')
# View(listEnrichrDbs())

dbs <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
         "CellMarker_2024", "Tabula_Muris",
         "Reactome_Pathways_2024")

# bg <- pvalsDf$gene

for(clust in seq(K)){
  message(clust)
  genes <- filter(sig_genes, cluster == clust) %>% pull(gene)
  enr <- enrichr(genes, dbs)
  for(db in dbs){
    if(nrow(enr[[db]])){
      message(db)  
      pdf(glue('output/enrichr/pam_shape_cluster/cluster_{clust}_{db}.pdf'), width=8, height = 4)
      plot(plotEnrich(enr[[db]], showTerms = 5, numChar = 50, 
                      # y_curves = "Count", 
                      # orderBy = "P.value", 
                      title = glue('Cluster-{clust} {db}'))) + 
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12)) # Change 12 to your desired font size
      
      dev.off()
    }
  }
  # printEnrich(enr, outFile = "txt", prefix = clust)
}

#--------------------------------
## make heatmap w/ DESeq tissue-level FC values
#--------------------------------
pdf(glue('figs/{CONTRAST}_LFC_heatmap_clustered_{K}.pdf'), width = 8, height = 10)
Heatmap(logFC[sig_genes$gene, ], right_annotation = row_anno,
        cluster_columns = T, cluster_rows = T,
        row_split = cluster_assignments,
        cluster_row_slices = F,
        name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=8))
dev.off()


myc <- colorRamp2(breaks = c(0, -log(.05), 5), colors = c("lightgray", "orange2", "darkred"))
pdf(glue('figs/{CONTRAST}_neglogp_heatmap_clustered_{K}.pdf'), width = 8, height = 10)
Heatmap(-log(p.vals)[sig_genes$gene, ], right_annotation = row_anno,
        cluster_columns = T, cluster_rows = T,
        col = myc,
        row_split = cluster_assignments,
        cluster_row_slices = F,
        name = 'neglogp',
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=8))
dev.off()
