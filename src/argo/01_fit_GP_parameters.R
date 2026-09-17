#' **01_fit_GP_parameters** : use vst-transformed bin-level pseudobulk counts data to fit GPs for each gene
#' we allow for age-specific noise terms, but assume common covariance structure that is learned jointly across all conditions
#' fitting a zero-mean GP for each gene, so we assume input data is centered
#' preserve original observation units by not scaling -- observed deviations (x1 - x2) indicate LFC values for the tested contrast, not linear count-deltas
#' 
#' @params KERNEL : allows fitting different covariance kernels to genes. available options are RBF, Matern3 (for matern 3/2 kernel) and Matern5 (for matern 5/2 kernel)
#' @params FIT_DONOR_SWITCH : allows flexible fitting of donor-level random intercepts
#' 
#' @outputs : kernelParams : fits the GPs for each gene and saves fitted parameters as a list
#' 

## run these to get the input data
source('src/scripts/00-load-packages.R')
# source('src/scripts/0B-load-data-for-modeling.R')
source('src/R/nll_and_gradients.R')

## Define the kernel to fit
KERNEL <- 'RBF'
## whether to fit donor-specific intercepts
FIT_DONOR_SWITCH = FALSE

## define cov kernel exponents
## pairwise distances and dist-squared
r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2

#--------------------------------
## set constraints for kernel bandwidth
#--------------------------------
min_bw <- min(dist(x_observed)) * 0.5
min_log_bw <- log(min_bw)
max_bw <- 3*diff(range(x_observed))
max_log_bw <- log(max_bw)

## list to store learned GP parameters for each gene type
kernelParams <- list()
for(gene in genes){
  #------------------------
  ## initialize parameter values
  #------------------------
  message(gene)
  ## abundance variable to be modeled
  y <- bin_meta %>% pull(gene)
  ## draw bandwidth parameter
  log_bw <- runif(1, log(.2), log(.5))
  
  ## initialize noise for donor effects
  ## I want this variable to capture the variance between donor mice at every point along the isodepth
  ## need to compute variance by mice in each bin and summarize
  ## however considering mice across age groups will confound with treatment effects -- 
  ## only compute within each age group and take mean
  log_sigma_d <- bin_meta %>% group_by(bin, age) %>%
    summarize(sd_val = sd(!!sym(gene), na.rm = T), .groups = "drop") %>%
    group_by(age) %>% # take median within age across all bins
    summarize(sd_age = mean(sd_val)) %>% # take the mean age-wise deviation for the gene
    pull(sd_age) %>% mean() %>% log()
  
  
  # compute empirical variation in gene by age group
  emp_sd <- bin_meta %>% group_by(age) %>% summarize(sd_val = sd(!!sym(gene))) %>% pull(sd_val)
  log_sigma_f <-  log(mean(emp_sd) * sqrt(0.7)) ## 70% variance to signal
  log_sigma_n <- log(emp_sd * sqrt(0.3)) ## 30% variance to noise
  parInt <- c(log_bw, log_sigma_f, log_sigma_d, log_sigma_n)
  
  #--------------------------------
  ## fit parameters
  #--------------------------------
  optimResult = stats::optim(par = parInt, 
                             fn = computeNegLogLikelihood, 
                             gr = computeNegLogLikelihoodGradient,
                             x = x_observed,
                             y = y,
                             P_obs = P_obs,
                             N = N, 
                             r = r_x,
                             r_squared = dsqx,
                             noise_diag = noise_diag,
                             slides_per_group = slides_per_age,
                             B = B,
                             cov_donor = cov_donor,
                             lambda = 1, # regularization for noise nuggets
                             kernel = KERNEL, # fit specific covariance-kernel to data
                             fit_donor = FIT_DONOR_SWITCH, # whether to decompose variance with donor-intercept
                             method = 'L-BFGS-B', 
                             lower = c(min_log_bw, rep(-3, 5)),
                             upper = c(max_log_bw,  rep(3, 5)),
                             control = list(maxit = 1000, trace = 2))
  
  estimatedParams = list(bandwidth = exp(optimResult$par[1]), 
                         sigma_f = exp(optimResult$par[2]), 
                         sigma_d = exp(optimResult$par[3]), 
                         sigma_n = exp(optimResult$par[4:6]))
  
  ## extract and store the parameters 
  kernelParams[[gene]] <- unlist(estimatedParams)
  
}

write_rds(kernelParams, glue('output/{KERNEL}_kernel_params.rds'))

#--------------------------------
## check fitted parameters
#--------------------------------
params <- do.call('rbind', kernelParams) 
colnames(params) <- c('bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')

pdf(glue('figs/{KERNEL}_fitted_params.pdf'), width = 6, height = 12)
Heatmap(log(params[, -3]), name = 'fitted GP \nparams', 
        row_names_gp = gpar(fontsize=4))
dev.off()

