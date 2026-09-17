#' plot genes in space that have large margnitude effects by funGP but are not includeded in any of the benchmarking genesets
#' list of genes are coming from the script *genesets-for-benchmarking.R*
#' I want to see if the effects are real. They seem to be from the posterior fits, but see them in space too


#===================================================================
## genes of interest that I can visualize
#===================================================================
library(Seurat)

CONTRAST <- 'old_vs_adult'
KERNEL <- 'RBF'

source('src/R/plot_conf_bands.R')
source('src/R/compute_exceedence_metrics.R')
source('src/R/compute_exceedence_statistics.R')

funGP_fits <- read_rds(glue('output/{KERNEL}_funGP_testing_{CONTRAST}_res.rds'))
load('data/temp/domain.rdata')

filtered.res <- res.fungp %>% arrange(desc(neglogp)) %>% filter(gene %in% genes & abs(deviation) > 1) 

stat_lfc <- compute_stat_continuous(funGP_fits, M)
pdf(glue('figs/FC_heatmap_select.pdf'), width = 8, height =6)
Heatmap(stat_lfc[filtered.res$gene, ], cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=12),
        column_names_gp = gpar(fontsize=4))
dev.off()


#===================================================================
## read list of seurat objects with representative images
#===================================================================
dxlist <- read_rds('~/Desktop/slide_list.rds')
dxlist <- lapply(dxlist, NormalizeData)
dxlist

## prepare common metadata file
dfm <- do.call(rbind, lapply(dxlist[3:6], function(dx) dx@meta.data))
dfm$donor_id <- droplevels(dfm$donor_id)
dfm$age <- droplevels(dfm$age)

#===================================================================
## load fitted kernel params
#===================================================================
kernelParams <- read_rds(glue('output/{KERNEL}_kernel_params.rds'))
params <- do.call('rbind', kernelParams)
colnames(params) <- c('bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')
FIT_DONOR_SWITCH = FALSE
## define cov kernel exponents
## pairwise distances and dist-squared
r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2

#===================================================================
## make posterior fits for genes i am visualizing
#===================================================================
# genes_to_plot <- c('Ctps1','Fezf2' ,'Myt1l', 'Cd19')
genes_to_plot <- filtered.res %>% arrange(desc(deviation)) %>% pull(gene)
condition <- 'diff_benchmark_multiple'
length(genes_to_plot)
source('src/scripts/02_visualize_fitted_GPs.R')

plots <- visualize_fitted_GP(genes_to_plot, params = params, FILTER_TO_CONTRAST = T)
px <- wrap_plots(plots, ncol = 4) + 
  plot_layout(guides = "collect")

ggsave(px,filename = file.path('figs', glue('posterior_fit_{condition}.pdf')),
       width=25,height=15)

#===================================================================
## visualize some representative genes in space
#===================================================================
## get roi for each gene
roi <- lapply(funGP_fits[genes_to_plot], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)

## extract counts for genes of interest across all slides
cts <- do.call('cbind', lapply(dxlist[3:6], function(dx) LayerData(dx)[genes_to_plot, ]))
isodepth_bins <- do.call('c', lapply(dxlist[3:6], function(dx) dx[['binned_isodepth', drop=T]]))


par(mai = c(1.0, 1, 0.8, 0.4)) 
boxplot(cts['Fezf2', ]~ dfm$age, xlab = 'Age group', ylab = 'Fezf2',
        font.axis = 2, font.lab = 2, cex.axis=1.5, cex.lab = 1.5)


## need to set thresholds across all three age groups by roi
highlight_roi <- function(dx, rx, gene_thresh, gene){
  expr <- LayerData(dx)[gene, ]
  alpha_map <- ifelse(dx$binned_isodepth %in% rx & expr, .5, 5e-3)
  size_map <-  ifelse(dx$binned_isodepth %in% rx, .5, 1e-2)
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
  gene_thresh <- quantile(cx[cx > 0], c(.25,  .75, 0.9));
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

pdf(glue('figs/{condition}_spatial_ROI_low_thresh.pdf'), width = 14, height = 15)
print(wrap_plots(plots, nrow=4) + 
        plot_annotation(title = glue('{CONTRAST}')))
dev.off()

