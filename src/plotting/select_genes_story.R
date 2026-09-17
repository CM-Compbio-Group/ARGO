library(ComplexHeatmap)
library(circlize)
library(Seurat)

setwd('/Volumes/umms-congma1/users/vravik/funGP-mouse_ageing/expression-modeling/')

source('src/scripts/00-load-packages.R')
source('src/scripts/0B-load-data.R')

KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'
source('src/R/plot_conf_bands.R')
source('src/R/compute_exceedence_metrics.R')
source('src/R/compute_exceedence_statistics.R')

kernel_params <- read_rds(glue('output/kernel_params_{KERNEL}.rds'))
kernel_params <- do.call('rbind', kernel_params)
colnames(kernel_params) <- c('bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')

## read funGP contrast results
exact_pvals <- read_rds(glue('output/funGP_testing_{KERNEL}_{CONTRAST}_res.rds'))
x_test <- read_rds('output/test_points.rds'); M <- length(x_test)

## compute exceedence statistics
exc.metrics <- do.call('rbind', lapply(exact_pvals, compute_exceedence_metrics, x_test))

## aggregate results from functional testing
pvalsDf <- exact_pvals %>%
  map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'gene') %>%
  bind_cols(exc.metrics) %>%
  mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         net_auec = auec_up - auec_down, 
         net_poe = poe_up + poe_down) %>%
  arrange(desc(stats)) %>% 
  mutate(Rank = row_number()) 

plot(pvalsDf$net_auec, pvalsDf$max_deviation)
plot(pvalsDf$Rank, pvalsDf$deviation)
range(pvalsDf$max_deviation[pvalsDf$max_deviation > 0])

#-------------------------------------
## visualize contrast test results
#-------------------------------------
stat_binary <- compute_stat_binary(exact_pvals, M)
stat_lfc <- compute_stat_continuous(exact_pvals, M)


filter(pvalsDf, max_deviation > .5 &  neglogp > -log(.05)) %>% 
  with(table(net_auec > 0))

goi <- filter(pvalsDf, max_deviation > .5 &  neglogp > -log(.05)) %>% pull(gene)

mx <- 1
myc <- colorRamp2(
  breaks = c(-mx, 0, mx), 
  colors = c("#377eb8", "gray60", "#e41a1c")
)


pdf(glue('figs/{CONTRAST}_FCmap_thresholded.pdf'), width = 8, height = 8)
Heatmap(stat_lfc[df$gene, ], cluster_columns = F, name = 'FC \nstatistic',
        show_row_names = T, col = myc,
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=4),
)
dev.off()

#-------------------------------------
#-------------------------------------
## no. of genes satisfying deviation based thresholds
thresh <- c(.15, .2, .25, .5, .75, 1)
map_dfr(thresh, function(t) {
  pvalsDf %>% filter(neglogp > -log(.05) & max_deviation > t & net_poe < 0.5) %>%
    summarise(auec_threshold = t, 
              upregulated = sum(net_auec > 0),
              downregulated = sum(net_auec < 0),
              total = n())
})


## filter to genes of interest
df <- filter(pvalsDf, net_poe < .5 & abs(deviation) > .5 &  neglogp > -log(.05))


pdf(glue('figs/{CONTRAST}_binary_heatmap_select.pdf'), width = 8, height = 6)
binary_colors <- c("1" = "red3","0" = "lightgrey", "-1" = "navy")
Heatmap(stat_binary[df$gene, ], cluster_columns = F, name = 'binary \nstatistic',
        show_row_names = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=4),
        col = binary_colors,
        heatmap_legend_param = list(color_bar = "discrete")
)
dev.off()

pdf(glue('figs/{CONTRAST}_FC_heatmap_select.pdf'), width = 8, height = 6)
Heatmap(stat_lfc[df$gene, ], cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=4),
        column_names_gp = gpar(fontsize=4))
dev.off()

#-------------------------------------
## pick genes to highlight
#-------------------------------------
# genes_to_plot <- c('Cdkn2a', 'Nfkbia', 'Apoe')
# genes_to_plot <- c('Satb2', 'Bcl11b', 'Tshz2', 'Cpne7', 'Sox10')
# genes_to_plot <- c('Zfp385b', 'Onecut2', 'Hmgcr', 'Pde10a', 'Igf2', 'Kcnip2', 'Olig1', 'Olig2', 'Sirt2')
# genes_to_plot <- c('Apoe', 'Abca1', 'Cxcl10', 'Hif3a', 'Crh', 'Slc30a3', 'Lrrc9', 'Hpse')
# genes_to_plot <- c('Adora2a',  'Hs3st4', 'Itga8', 'Clec10a', 'Sulf1', 'Cldn34c1')
# genes_to_plot <- c('Serpina3n', 'C1qa', 'Tgfb1', 'Itgb2', 'Foxj1' )
# genes_to_plot <- c('Hs3st4', 'Npr3', 'Pvalb', 'Opalin', 'Eomes', 'Igf2', 'C1qa', 'Gsdmd', 'Car3', 'Scube1')

genes_to_plot <- c('Serpina3n', 'Lmo1', 'Adora2a', 'Hs3st4')

genes_to_plot <- rownames(mask_discordant)[1:5]
condition <- 'discordant_trends'
length(genes_to_plot)

source('src/scripts/02_visualize_fitted_GPs.R')
px <- wrap_plots(plots, ncol = 5) + 
  plot_layout(guides = "collect")

ggsave(px,filename = file.path('figs', glue('{CONTRAST}_posterior_fit_{condition}.pdf')),
       width=22,height=8)

conf_band_plots <- lapply(genes_to_plot, function(gene) 
  plot_CI_bands(gene, x = x_test, 
                band = exact_pvals[[gene]][['band']],
                stat = exact_pvals[[gene]][['statistic']], 
                alpha = 95)
)

p <- wrap_plots(conf_band_plots, ncol = 5) + 
  plot_layout(guides = "collect")

ggsave(p,filename = file.path('figs', glue('{CONTRAST}_CI_bands_{condition}.pdf')),
       width=22,height=8)

#-----------------------------------------------
## make spatial plots with two genes per type
#-----------------------------------------------
## specify genes to view
genes_to_plot <- c('Serpina3n', 'Tgfb1','Foxj1', #white-matter specific
             'Apoe', 'Abca1','Slc30a3', # cortical layer speciifc
             'Adora2a','Sulf1', 'Itga8', # striatum specific
             'Igf2', 'Olig1', 'Olig2', #sub-ventricular region depletion
             'Satb2', 'Bcl11b', 'Tshz2' #cortical senescence markers
)

## read list of seurat objects with representative images
dxlist <- read_rds('~/Desktop/slide_list.rds')
dxlist <- lapply(dxlist, NormalizeData)

## align the two images
dx<-dxlist[[2]]
tc <- dx@reductions$X_spatial_coords@cell.embeddings
tc[,2] <- -tc[,2]
dx@reductions$X_spatial_coords@cell.embeddings <- tc
DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'tissue')
dxlist[[2]] <- dx

#-----------------------------------------------
## get their domain maps
#-----------------------------------------------
figs <- lapply(dxlist, function(dx) {
  DimPlot(dx, reduction = 'X_spatial_coords', group.by = 'tissue') + 
    ggtitle('') + 
    labs(x = '', y = '') + 
    theme_minimal() + 
    theme( panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
      axis.text = element_blank()) 
})

p <- wrap_plots(figs, ncol = 2) + 
  plot_layout(guides = "collect")

ggsave(p,filename = file.path('figs', glue('{CONTRAST}_spatial_domains.pdf')),
       width=12,height=4)

#-----------------------------------------------
#-----------------------------------------------
## get roi for each gene
roi <- lapply(exact_pvals[genes_to_plot], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)

## extract counts for genes of interest across all slides
cts <- do.call('cbind', lapply(dxlist[3:6], function(dx) LayerData(dx)[genes_to_plot, ]))
isodepth_bins <- do.call('c', lapply(dxlist[3:6], function(dx) dx[['binned_isodepth', drop=T]]))
domainx <- dfm[colnames(cts), ] %>% pull(tissue)

## need to set thresholds across all three age groups by roi
highlight_roi <- function(dx, rx, gene_thresh, gene){
  expr <- LayerData(dx)[gene, ]
  alpha_map <- ifelse(dx$binned_isodepth %in% rx & expr, .25, 5e-3)
  size_map <-  ifelse(dx$binned_isodepth %in% rx, .25, 1e-2)
  p <- FeaturePlot(dx, feature = gene,
                   reduction = 'X_spatial_coords', 
                   # min.cutoff = gene_thresh[1],
                   # max.cutoff = gene_thresh[2],
                   alpha = alpha_map,
                   pt.size=size_map) + 
    scale_color_gradientn(
      colors = c("lightgray",  "red2", "darkred"),
      limits = gene_thresh[c(1,3)],
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
    gene_thresh <- quantile(cx[cx > 0], c(.5,  .75, 0.99));
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

pdf(glue('figs/{condition}_spatial_ROI.pdf'), width = 12, height = 30)
print(wrap_plots(plots, nrow=10) + 
        plot_annotation(title = glue('{CONTRAST}')))
dev.off()


# for(i in seq(5)){
#   index = 3*i - 2
#   message(index)
#   pdf(glue('figs/{CONTRAST}_select_genes_spatial_{i}.pdf'), width = 10, height = 10)
#   p <- wrap_plots(figs[index : (3*i)], nrow = 3) + 
#     plot_annotation(title = glue('{CONTRAST}'))
#   print(p)
#   dev.off()
# }



#-----------------------------------------------
## how many of these genes are found uniquely by fungp?
#-----------------------------------------------
df2 <- vroom(glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv')) %>% arrange(desc(stats))
df2 <- df2 %>% mutate(neglogp = pmax(neglogp_combined, neglogp_min)) %>% 
  filter(neglogp > -log(.05) & abs(maxlogFC) > 0.25)

goi <- setdiff(df$gene, df2$gene)
pdf(glue('figs/{CONTRAST}_FC_heatmap_fungp_unique.pdf'), width = 8, height = 6)
Heatmap(stat_lfc[goi, ], cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=12),
        column_names_gp = gpar(fontsize=4))
dev.off()

filter(df, gene %in% goi)

#------------------------------------
## find discordant genes in top ranks
#------------------------------------
## ordering of deseq results depends on pval statistic we use
plot(rank(df2$neglogp_combined), rank(df2$neglogp_min))

# list of discordant genes will vary
goi <- setdiff(pvalsDf %>% pull(gene) %>% head(n = 50), 
        df2 %>% pull(gene) %>% head(n = 50))

filter(pvalsDf, gene %in% goi)

