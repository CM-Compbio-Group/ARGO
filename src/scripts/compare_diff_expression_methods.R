## gather all results and make plots
setwd('/Volumes/umms-congma1/users/vravik/funGP-mouse_ageing/expression-modeling/')
source('src/scripts/00-load-packages.R')

KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'

res.deseq <- read_rds('output/deseq_DGE_res.rds')[['dge']]
res.edger <- read_rds('output/edgeR_DGE_res.rds')[['dge']]
res.limma <- read_rds('output/limma_DGE_res.rds')[['dge']]
res.tradeseq <- read_rds('output/tradeSeq_DGE_res.rds')
res.fungp <- vroom(glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv')) 
res.lamian <- read_rds('output/lamian_DGE_res.rds')

res.fungp %>%  ggplot(aes(x = deviation, y = neglogp)) + geom_point()

## extract statistics
res.lamian <- res.lamian$statistics %>% 
  rownames_to_column('gene') %>% 
  as_tibble() %>%
  mutate(neglogp = -log.pval.overall) %>% 
  arrange(desc(neglogp)) %>% mutate(Rank = row_number())


res.tradeseq <- res.tradeseq %>% arrange(desc(waldStat_conds2vs3)) %>% 
  mutate(Rank = row_number(), 
        neglogp = -log(pvalue_conds2vs3 + 1e-16)) %>% as_tibble()

res.list <- list('deseq' = res.deseq,
                 'edgeR' = res.edger,
                 'limma' = res.limma,
                 'tradeseq' = res.tradeseq,
                 'lamian' = res.lamian,
                 'funGP' = res.fungp)

res.list

write_rds(res.list, glue('output/{CONTRAST}_DGE_res_list.rds'))

#====================================================#
## run enrichment analysis against genesets
#====================================================#
## filter to thresholded genes
pval_thresh <- 0.01
res.thresh <- lapply(res.list, function(df) filter(df, neglogp > -log(pval_thresh)))
sapply(res.thresh, nrow)

library(enrichR)
dbs <- 'GTEx_Aging_Signatures_2021'
dbs <- c("GO_Biological_Process_2026", "TISSUES_Curated_2025", "MSigDB_Hallmark_2020",
         "CellMarker_2024", "GTEx_Aging_Signatures_2021", "Tabula_Muris",
         "Reactome_Pathways_2024")

dir.create(glue('output/enrichr/diffx/{CONTRAST}_comparative_w_bg'))

for(method in names(res.thresh)){
  message(method)
  genes <- res.thresh[[method]]$gene
  enr <- enrichr(genes, dbs, background = res.fungp$gene)
  
  for(db in dbs){
    if(nrow(enr[[db]])){
      message(db)  
      pdf(glue('output/enrichr/diffx/{CONTRAST}_comparative_w_bg/{method}_{db}.pdf'), width=10, height = 6)
      plot(plotEnrich(enr[[db]], showTerms = 10, numChar = 50, 
                      y  = "Count", 
                      title = glue('{method} {db}'))) + 
        theme(axis.text.y = element_text(size = 16)) # Change 12 to your desired font size
      
      dev.off()
    }
  }
}

## what are the aging-up genes in lamian?
head(enr[[db]])

genes <- 'C4B;IFITM3;FEZF2;GFAP;FAM107A'
genes <- str_split(genes, ';')[[1]] %>% str_to_title()
genes
res.fungp %>% filter(gene %in% genes)

#====================================================#
## compute the auprc curves
#====================================================#
#-------------------------------------------------------------------------
## aging-signature genes from https://www.cell.com/cell/fulltext/S0092-8674%2823%2900805-X
#-------------------------------------------------------------------------
cas <- vroom('data/sheets/cas-genes.txt', delim = '\t')

p1 <- get_ranked_dge_plot(res.fungp, method = 'funGP', markers = cas$`Gene Symbol`)
p2 <- get_ranked_dge_plot(res.lamian, method = 'lamian', markers = cas$`Gene Symbol`)
pdf('figs/CAS-rankes-statistics.pdf', width = 12)
wrap_plots(p1 | p2)
dev.off()

## find direction of change for these genes in all tests
# but i am not keeping track of deviation anymore

filter(res.lamian, gene %in% cas$`Gene Symbol`)
filter(res.fungp, gene %in% cas$`Gene Symbol`)

cas <- filter(res.fungp, gene %in% cas$`Gene Symbol` & deviation > 0) %>% pull(gene)

#----------------------------------------------------------------------------------------
## list of age-related differential genes from https://www.nature.com/articles/s41586-024-08350-8#Sec8
#----------------------------------------------------------------------------------------
library(readxl)
library(PRROC)
## filtered list of genes with only restricted celltypes
neuronal.dge <- read_xlsx('data/sheets/cortex_striatum_age_DE_filtered.xlsx', sheet = 1)
glial.dge <- read_xlsx('data/sheets/cortex_striatum_age_DE_filtered.xlsx', sheet = 2)
names(neuronal.dge)[7] <- names(glial.dge)[7] <- 'label'
age.dge <- rbind(neuronal.dge, glial.dge)

# restrict to MERFISH panels
age.dge <- filter(age.dge, gene %in% res.fungp$gene)
## filter out genes with dual trends across cell types
gene_summary <- age.dge %>% group_by(gene) %>%
  summarize(upreg = sum(age_effect_size > 0),
            downreg = sum(age_effect_size < 0),
            mixed = any(age_effect_size > 0) & any(age_effect_size < 0)
  ) %>% arrange(desc(mixed))

sum(gene_summary$mixed)
sum(gene_summary$upreg > 1)
sum(gene_summary$downreg > 1)

filter(gene_summary, gene %in% CAS$`Gene Symbol`)

## can filter out the mixed-signal genes if necessary
## can threshold by how many cells have differential signature
n_thresh <- 1

# markers <- filter(gene_summary, upreg > n_thresh | downreg > n_thresh) %>% pull(gene)


#------------------------------------------------------
## make the auprc curves
#------------------------------------------------------
compute_ranked_auprc <- function(ranked_list, signature){
  binary_labels <- ranked_list %in% signature
  scores <- (seq(length(ranked_list), 1))
  pg_score <- scores[binary_labels]
  ng_score <- scores[!binary_labels]
  pr_obj <- pr.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)
  return(pr_obj)
}


# markers <- gene_summary %>% pull(gene)
markers <- cas
# signature <- CAS$`Gene Symbol`
## Rank is based on directional statistics!
## I only want to rank by p-values!!
pr.curves <- lapply(res.list, function(df)  {
      compute_ranked_auprc(df %>% arrange(desc(neglogp)) %>% pull(gene), signature = markers)
    }
  )

# SX <- 'Jin et al aging markers'
SX <- 'Common Aging Signature'
cols <- c( 'red3', 'lightskyblue2', 'lightgreen', 'lightpink2', 'lightsalmon', 'orchid2')

pdf(glue('figs/{CONTRAST}_AUPRC_comparative_{SX}_filtered.pdf'), width = 8, height = 8)
par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
plot(pr.curves[['funGP']], col = adjustcolor(cols[1], alpha.f = 1), main = glue("{SX} \n Precision-Recall Curve"), 
     auc.main = F, lty=2)
## plot other curves
for(ii in seq(1,5)){
  method = names(pr.curves)[ii]
  plot(pr.curves[[method]], col = adjustcolor(cols[ii + 1], alpha.f = .5), lty = 2, lwd = 2, add = TRUE)
}
## make legends
lgx <- paste("FunGP (AUPRC = ", round(pr.curves[['funGP']]$auc.integral, 3), ")")
for(method in names(pr.curves)[1:5]){
  lgx <- c(lgx, 
              glue("{method} (AUPRC = {round(pr.curves[[method]]$auc.integral, 3)}) ")
              )
}
legend("topright", cex=.8,
       legend = lgx,
       col = cols, 
       lty = 2,
       lwd = 2,             # Line width
       bty = "n")           # Removes the box border around the legend (optional)

dev.off()

#====================================================#
## find genes uniquely identified by funGP
#====================================================#
bin_dds <- read_rds('data/temp/bin_dds.rds')
setdiff(res.thresh$lamian$gene, res.thresh$funGP$gene)
plotCounts(bin_dds, 'C1qa', intgroup = 'age_bin')

diff_genes <- lapply(res.thresh, '[[', 'gene')
## set of genes uniquely identified by funGP
dgx <- lapply(diff_genes[1:5], function(z) setdiff(diff_genes[['funGP']], z))
sapply(dgx, length)

gx1 <- Reduce(intersect, dgx)
gx2 <- Reduce(intersect, dgx[-4])

## order genes by their magnitude
res.fungp %>% filter(gene %in% gx1)  %>% ggplot(aes(x = deviation, y = neglogp)) + geom_point()
goi <- res.fungp %>% filter(gene %in% gx1) %>% {bind_rows(slice_head(., n = 3), slice_tail(., n = 3))}

res.fungp %>% filter(gene %in% gx2)  %>% ggplot(aes(x = deviation, y = neglogp)) + geom_point()
goi <- res.fungp %>% filter(gene %in% gx2) %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))}

  
## visualize trends for these genes
my_colors <- colorRamp2(
  breaks = c(-1, 0, 1), 
  colors = c("blue", "white", "red")
)
Heatmap(stat_lfc[goi$gene, ], cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=12), col = my_colors,
        column_names_gp = gpar(fontsize=4))

goi.res <- lapply(res.list, function(df) filter(df, gene %in% goi) %>% arrange(match(gene, goi)))
#-------------------------------------
## visualize genes in space
#-------------------------------------
## assign regions from discordant mask
library(Seurat)
genes_to_plot <- goi$gene
condition <- 'uniq_high_deviation'
length(genes_to_plot)

source('src/scripts/02_visualize_fitted_GPs.R')
px <- wrap_plots(plots, ncol = 3) + 
  plot_layout(guides = "collect")

ggsave(px,filename = file.path('figs', glue('{CONTRAST}_posterior_fit_{condition}.pdf')),
       width=15,height=8)
#-----------------------------------------------
#-----------------------------------------------
conf_band_plots <- lapply(genes_to_plot, function(gene) 
  plot_CI_bands(gene, x = x_test, 
                band = exact_pvals[[gene]][['band']],
                stat = exact_pvals[[gene]][['statistic']], 
                alpha = 95)
)

p <- wrap_plots(conf_band_plots, ncol = 3) + 
  plot_layout(guides = "collect")

ggsave(p,filename = file.path('figs', glue('{CONTRAST}_CI_bands_{condition}.pdf')),
       width=15,height=8)

#-----------------------------------------------
#-----------------------------------------------
dxlist <- read_rds('~/Desktop/slide_list.rds')
dxlist <- lapply(dxlist, NormalizeData)


## make spatial plots for tissue boundaries and isodepth
figs <- list()
for(dx in dxlist[4:5]){
  donor <- as.character(dx$donor_id[1])
  figs[[donor]] <- FeaturePlot(dx, reduction = 'X_spatial_coords', features  = 'isodepth',
                               min.cutoff = 0, max.cutoff = 1) +
    scale_color_gradient(low = 'lightgray',high = 'darkred', breaks = seq(0,1,.2)) +
  # figs[[donor]] <- DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'tissue') + 
    ggtitle('') +
    labs(x = '', y = '') + 
    theme_minimal(base_family = 'Arial') + 
    theme(axis.text = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

image_plots <- wrap_plots(figs, ncol = 2) + 
  plot_layout(guides = "collect") 

pdf(glue('figs/isodepth_select.pdf'), width = 10, height = 4)
print(image_plots)
dev.off()



## find rois for the genes
roi <- lapply(exact_pvals[genes_to_plot], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
(lx <- sapply(roi, length))

## extract metadata from all slides
dfm <- bind_rows(lapply(dxlist, function(dx) dx@meta.data))
dfm <- filter(dfm, age %in% c('24wk', '90wk'))

## extract counts for genes of interest across all slides
cts <- do.call('cbind', lapply(dxlist[3:6], function(dx) LayerData(dx)[genes_to_plot, ]))
isodepth_bins <- do.call('c', lapply(dxlist[3:6], function(dx) dx[['binned_isodepth', drop=T]]))

## need to set thresholds across all three age groups by roi
highlight_roi <- function(dx, rx, gene_thresh, gene){
  expr <- LayerData(dx)[gene, ]
  alpha_map <- ifelse(dx$binned_isodepth %in% rx & expr, .25, 5e-3)
  size_map <-  ifelse(dx$binned_isodepth %in% rx, .25, 1e-2)
  p <- FeaturePlot(dx, feature = gene,
                   reduction = 'X_spatial_coords', 
                   alpha = alpha_map,
                   pt.size=size_map) + 
    scale_color_gradientn(
      colors = c("lightgray",  "red3", "darkred"),
      limits =  range(gene_thresh),
      breaks = gene_thresh
    ) + 
    ggtitle('') + 
    labs(x = '', y = '') + 
    theme_minimal(base_family = 'Arial') + 
    theme(axis.text = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}

figs <- lapply(genes_to_plot, function(gene) {
  message(gene)
  rx <- roi[[gene]];
  cx <- cts[gene, (isodepth_bins %in% rx)];
  gene_thresh <- quantile(cx[cx > 0],  c(.5, .75, 0.95, .995));
  lapply(dxlist[3:6], highlight_roi, rx, gene_thresh, gene )
  }
)

names(figs) <- genes_to_plot
plots <- lapply(genes_to_plot, function(gene) {
  message(gene)
  row_plot <- wrap_plots(figs[[gene]], ncol = 4) + 
    plot_layout(guides = "collect") + 
    plot_annotation(title = glue('{gene} ROI'))
  wrap_elements(row_plot)
}
)

pdf(glue('figs/{condition}_select_spatial_ROI.pdf'), width = 12, height = 12)
print(wrap_plots(plots, nrow=length(genes_to_plot)) + 
        plot_annotation(title = glue('{CONTRAST}')))
dev.off()

#-------------------------------------
#-------------------------------------
## make domain plots from deseq with same pvalue thresholds
deseq <- read_rds('output/limma_DGE_res.rds')
p.vals <- deseq$p.vals[genes_to_plot, ]
p.vals[p.vals < pval_thresh]

