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

## for parallel processing
library(futurize)

## for fitting scDesign and simulating data
source('testing/src/R/fit_scDesign3.R')
source('testing/src/R/generate_perturbed_data.R')
source('testing/src/R/generate_interval_perturbation.R')

## differential testing
source('testing/src/R/run_differential_testing.R')
source('testing/src/R/run_interval_evaluation.R')
source('testing/src/R/utils.R')

## for baseline methods
source('testing/src/R/get_cct_neglogp.R')
source('src/R/summarize_dge_results.R')
# 
## for argo
source('testing/src/R/run_argo.R')
source('testing/src/R/simple_likelihood.R')
source('testing/src/R/simple_contrast_test_with_funGP.R')
source('testing/src/R/generate_cc_matrix.R')
source('testing/src/R/compute_analytical_pvalue.R')
source('testing/src/R/derive_KLT_CI_bands.R')
source('testing/src/R/compute_exceedence_statistics.R')




