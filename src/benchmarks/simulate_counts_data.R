#' generate single cell data for benchmarking argo
# devtools::install_github("SONGDONGYUAN1994/scDesign3")
library(scDesign3)
library(tidyverse)
library(SingleCellExperiment)
library(useful)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())


## I am going to use the pancreatic endocrinogenesis dataset as baseline
# download.file(url = "https://figshare.com/ndownloader/files/40581992", 
#               destfile = 'data/pancreas_sce.rds')

pancreas_sce <- read_rds('data/PANCREAS_sce.rds')
print(pancreas_sce)

## reduce to 100 genes
pancreas_sce <- pancreas_sce[1:200, ]
plot_reduceddim(ref_sce = pancreas_sce, 
                sce_list = list(),
                name_vec = c("Reference"),
                assay_use = "logcounts", 
                if_plot = TRUE, 
                color_by = "pseudotime", 
                n_pc = 20)

assayNames(pancreas_sce)
head(as.data.frame(colData(pancreas_sce)))

keep_genes <- Matrix::rowSums(counts(pancreas_sce) > 0) >= 20
pancreas_sce <- pancreas_sce[keep_genes, ]
plot(
  colData(pancreas_sce)$pseudotime,
  Matrix::colSums(counts(pancreas_sce)),
  xlab = "Pseudotime",
  ylab = "Library size"
)


set.seed(148)

PANCREAS_data <- construct_data(
  sce = pancreas_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = "pseudotime",
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "1"
)

PANCREAS_marginal <- fit_marginal(
  data = PANCREAS_data,
  predictor = "gene",
  mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
  sigma_formula = "1", 
  family_use = "nb",
  n_cores = 16,
  usebam = FALSE
)

PANCREAS_copula <- fit_copula(
  sce = pancreas_sce,
  assay_use = "counts",
  marginal_list = PANCREAS_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 16,
  input_data = PANCREAS_data$dat,
)


PANCREAS_para <- extract_para(
  sce = pancreas_sce,
  marginal_list = PANCREAS_marginal,
  n_cores = 12,
  family_use = "nb",
  new_covariate = NULL,
  data = PANCREAS_data$dat
)

#===================================================================
## inject condition-specific effects into the mean parameters
#===================================================================
mu_A <- PANCREAS_para$mean_mat
mu_B <- mu_A

t <- PANCREAS_data$dat$pseudotime

## sample different patterns of condition-specific effects across the genes
# Randomly assign 25 genes to each pattern
genes <- rownames(counts(pancreas_sce))
gene_order <- sample(genes)
G <- length(genes)


# 25 null genes; 75 genes with Gaussian-bump effects
gene_order <- sample(genes)
n_null <- G / 4 # quarter of the genes are unaffected
null_genes <- gene_order[1:n_null]
effect_genes <- gene_order[(n_null + 1):G]

lfc_mat <- matrix(
  0,
  nrow = nrow(mu_A),
  ncol = G,
  dimnames = dimnames(mu_A)
)

# Gene-specific effect parameters
effect_info <- data.frame(
  gene = effect_genes,
  amplitude = runif(length(effect_genes), 0.5, 1.5),
  center = runif(length(effect_genes), 0.2, 0.8),
  bandwidth = sample(
    c(0.025,0.05, 0.10, 0.20),
    length(effect_genes),
    replace = TRUE
  )
)

for (j in seq_len(nrow(effect_info))) {
  g <- effect_info$gene[j]
  
  lfc_mat[, g] <-
    effect_info$amplitude[j] *
    exp(
      -(t - effect_info$center[j])^2 /
        (2 * effect_info$bandwidth[j]^2)
    )
}

mu_B <- mu_A * exp(lfc_mat)

matplot(
  t,
  lfc_mat[, effect_info$gene[1:10]],
  lty = 2, lwd = .5,
  xlab = "Pseudotime",
  ylab = "True LFC"
)


#===================================================================
## simulate counts data
#===================================================================
## I need to simulate data for multiple samples from each condition
n_samples <- 3 # no. of observations per condition

# Condition A: unperturbed fitted means
para_A <- PANCREAS_para
# Condition B: perturbed fitted means
para_B <- PANCREAS_para
para_B$mean_mat <- mu_B

generate_counts <- function(para){
  simu_new(sce = pancreas_sce,
           mean_mat = para$mean_mat,
           sigma_mat = para$sigma_mat,
           zero_mat = para$zero_mat,
           quantile_mat = NULL,
           copula_list = PANCREAS_copula$copula_list,
           n_cores = 16,
           family_use = "nb",
           input_data = PANCREAS_data$dat,
           new_covariate = PANCREAS_data$newCovariate,
           important_feature = PANCREAS_copula$important_feature,
           filtered_gene = PANCREAS_data$filtered_gene)
}

counts_A <- vector("list", n_samples)
counts_B <- vector("list", n_samples)

for(r in seq_len(n_samples)){
  set.seed(100 + r)
  counts_A[[r]] <- generate_counts(para_A)
  counts_B[[r]] <- generate_counts(para_B)
}


truth <- list(
  lfc_mat = lfc_mat,
  effect_info = effect_info,
  null_genes = null_genes, 
  pseudotime = PANCREAS_data$dat$pseudotime
)


sim_counts <- do.call(cbind, c(counts_A, counts_B))
## append unique sample name for each datapoint
sim_meta <- data.frame(condition = factor(rep(c('condition_A', 'condition_B'), 
                                          each = n_samples * ncol(pancreas_sce))),
                       replicate = rep(factor(rep(paste('replicate', 1:3, sep = '_'), 
                                              each = ncol(pancreas_sce))), 
                                       2), 
                       sample = rep(factor(paste('sample', LETTERS[1:6], sep = '-')), 
                                    each = ncol(pancreas_sce)),
                       pseudotime = rep(PANCREAS_data$dat$pseudotime, 2*n_samples)) %>% as_tibble()


save(sim_counts, sim_meta, truth, file = 'data/scdesign_simulation.rdata')
