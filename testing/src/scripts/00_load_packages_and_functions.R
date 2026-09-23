## load necessary libraries for the experiment

library(scDesign3)
library(glue)
library(tidyverse)
library(SingleCellExperiment)
library(useful)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(pheatmap)
library(ComplexHeatmap)

library(DESeq2)
library(edgeR)
library(limma)
library(tradeSeq)
library(Lamian)
require(patchwork)
require(zeallot)

library(readxl)
library(PRROC)
library(extrafont)

# ## for baseline methods
# source('src/R/get_cct_neglogp.R')
# source('src/R/summarize_dge_results.R')
# source('src/R/get_ranked_dge_plot.R')
# 
# ## for argo
# source('src/R/simple-likelihood.R')
# source('src/R/simple_contrast_test_with_funGP.R')
# source('src/R/generate_cc_matrix.R')
# source('src/R/compute_analytical_pvalue.R')
# source('src/R/derive_KLT_CI_bands.R')
# source('src/R/plot_conf_bands.R')
# source('src/R/compute_exceedence_metrics.R')
# source('src/R/compute_exceedence_statistics.R')
# 
# source('src/benchmarks/run_baseline_testing.R')
# source('src/benchmarks/run_argo.R')
# source('src/benchmarks/compute_eval_metrics.R')
