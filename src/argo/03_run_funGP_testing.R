#' **03_run_funGP_testing**:  run contrast tests w/ learned GP fits
#' @params CONTRAST_VECTOR : defines the contrast to test
#' @params KERNEL : what kernel is used for likelihood computation
#' @params KLT_CONF_BAND : whether to fit confidence bands for the gene, using monte carlo samples
#' @params genes_to_test : list of genes to test. defaults to all genes in the panel 
#' @params n_mc : number of monte-carlo simulations to derive level-alpha confidence bands
#' saves list of funGP testing objects for every gene to rds object in file

#--------------------------------
## run differential testing with funGP
#--------------------------------
# source('src/scripts/00-load-packages.R')
source('src/R/generate_cc_matrix.R')
source('src/R/compute_analytical_pvalue.R')
source('src/R/derive_KLT_CI_bands.R')
source('src/R/contrast-test-with-funGP.R')

## load fitted kernel parameters
kernelParams <- read_rds(glue('output/{KERNEL}_kernel_params.rds'))
params <- do.call('rbind', kernelParams)
colnames(params) <- c('bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')

## control variables
KERNEL <- 'RBF'
FIT_DONOR_SWITCH <- FALSE
KLT_CONF_BAND <- TRUE
CONTRAST <- 'old_vs_adult'
CONTRAST_VECTOR <- c(0, -1, 1)
genes_to_test <- rownames(params)


funGP_fits <- list()
for(gene in genes_to_test){
  message(gene)
  y <- bin_meta %>% pull(gene)
  funGP_fits[[gene]]  <- contrast_test_with_funGP(contrast_vector = CONTRAST_VECTOR, 
                                                   y = y, par= params[gene, ], 
                                                   r_x = r_x, r_t = r_t, r_tx = r_tx,
                                                   dsqx = dsqx, dsqt = dsqt, dsqtx = dsqtx,
                                                   x_observed = x_observed,
                                                   x_test = x_test,
                                                   P_obs = P_obs,
                                                   W = W, B = B, 
                                                   N = N, M = M,
                                                   noise_diag = noise_diag, 
                                                   kernel = KERNEL, 
                                                   fit_donor = FIT_DONOR_SWITCH, 
                                                   cov_donor = cov_donor,
                                                   return_conf_band = KLT_CONF_BAND,
                                                   alpha = 0.95, 
                                                   n_mc = 1e4
  )
}

## save diff testing results to file
write_rds(funGP_fits, glue('output/{KERNEL}_funGP_testing_{CONTRAST}_res.rds'))

