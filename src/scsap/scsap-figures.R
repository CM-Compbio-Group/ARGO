#' benchmarking related plots for poster

library(SingleCellExperiment)
library(glue)
library(tidyverse)
library(scDesign3)
library(patchwork)

dataset <- 'PANCREAS'
sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))
# add pseudotime bins to object
sce$bin <- factor(ntile(sce$pseudotime, n = 50))

p1 <- plot_reduceddim(ref_sce = sce, 
                sce_list = list(),
                name_vec = 'Pseudotime',
                assay_use = "logcounts", 
                if_plot = TRUE, 
                color_by = c("pseudotime"), 
                n_pc = 20)[[2]]

p2 <- plot_reduceddim(ref_sce = sce, 
                      sce_list = list(),
                      name_vec = 'Cell Type',
                      assay_use = "logcounts", 
                      if_plot = TRUE, 
                      color_by = c("cell_type"), 
                      n_pc = 20)[[2]] + 
  labs(color = 'Cell Type') + 
  theme(legend.position = "right")

p3 <- plot_reduceddim(ref_sce = sce, 
                      sce_list = list(),
                      name_vec = 'Pseudotime Bins',
                      assay_use = "logcounts", 
                      if_plot = TRUE, 
                      color_by = c("bin"), 
                      n_pc = 20)[[2]] +  
  theme(legend.position = "none")




p <- p1 | p3 | p2
ggsave(plot = p, filename = 'figs/pancreas_UMAP.pdf', width = 12, height = 4)


#========================================================================
## visualize effects for the genes
#========================================================================
zeallot::`%<-%`(c(sce, sce.dat, sce.copula, sce.para), sce_fits) #unlist learned variables for dataset
effect_genes <- sample(rownames(sce), size = 100, replace = T)
unique_effect_genes <- make.unique(effect_genes)

t_ <- sce.dat$newCovariate$pseudotime
ord <- order(t_)

## extract the mean parameters for these genes
mu_A <- sce.para$mean_mat[, effect_genes]
lfc_mat <- matrix(
  0,
  nrow = nrow(mu_A),
  ncol = ncol(mu_A),
  dimnames = dimnames(mu_A)
)

# amplitude and bandwidth are fixed, only vary where the perturbation occurs, and sign 
effect_info <- tibble(
  gene = unique_effect_genes,
  amplitude = sample(c(0.3, 0.6, 1), size = length(effect_genes), replace = T), 
  bandwidth = sample(c(0.03, 0.06, 0.1), size = length(effect_genes), replace = T),
  sign = sample(c(-1, 1), size = length(effect_genes), replace = T),
  center = runif(length(effect_genes), 0.2, 0.8)
)

print(table(effect_info$amplitude, effect_info$bandwidth))

## inject perturbations on mean parameter
for (j in seq_along(effect_genes)) {
  lfc_mat[, j] <-
    effect_info$sign[j] * effect_info$amplitude[j] *
    exp(
      -(t_ - effect_info$center[j])^2 /
        (2 * effect_info$bandwidth[j]^2)
    )
}

pdf('figs/representative_perturbations.pdf', width = 8, height = 6)
par(mai = c(1.2,1.2,1.2,.8))
matplot(t_[ord], lfc_mat[ord, 1:10], type = 'l', lwd = 2, cex.axis = 1.25, cex.lab = 1.25, cex.main =1.5,
        xlab = 'Pseudotime', ylab = 'LFC values', main = 'Representative perturbations between conditions')
dev.off()
