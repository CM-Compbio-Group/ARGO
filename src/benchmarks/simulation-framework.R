#' simulation framework for testing how good our algorithm is in detecting differential expt1 expression trends
#' simulate 2-condition scenario. $H_0 : u_1 = u_2, H_1 : \exists x \in X s.t. u_1(x) != u_2(x) $
#' execute this by drawing two latent functions $u_1, u_2 \in GP(0, \Sigma)$, and some time domain where we add introduce the deviation
#' make baseline curve, perturb with continuous function over interval for alternate condition
#' make the simulated datasets in this script

#' set parameter values empirically from observations
#' simulate data over defined combination of perturbation strength and interval extents

# setwd('/Volumes/umms-congma1/users/vravik/funGP-mouse_ageing/ARGO/')
source('src/scripts/00-load-packages.R')
library(MASS)
#---------------------------------------------------------------
## set up simulation framework
## set dimensions of domain
# draw noise and bandwidth parameters to get a new baseline curve
## define parameters of covariance kernel
## have a hierarchical model to pick kernel parameters
#---------------------------------------------------------------
KERNEL <- 'RBF' # cov kernel we are using for generating the latent processes

## load fitted parameter values for mouse brain aging data to get kernel and noise-parameter estimates
fittedParams <- read_rds(glue('output/{KERNEL}_kernel_params.rds'))
fittedParams <- do.call('rbind', fittedParams) %>% data.frame() %>% rownames_to_column('gene') %>% as_tibble()
colnames(fittedParams) <- c('gene','bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')
fittedParams

## get ranges for noise parameters
param_range <- sapply(fittedParams %>% dplyr::select(bandwidth, sigma_f, sigma_n2, sigma_n3), quantile, probs = c(.25, .75))
#---------------------------------------------------------------
#---------------------------------------------------------------
T <- 100
domain_x <- seq(0, 1, length.out = T)
r_squared <- as.matrix(dist(domain_x)^2)
N_CURVES <- 10

## set combination of values for perturbations
perturb_strength_list <- c(.05, .1, .2, .3, .5, 1, 2) # as a percent of cov-kernel noise
perturb_interval_list <- c(.05, .1, .2, .3,.4,.5) # extent of interval in domain
(N_EXPERIMENTS <- length(perturb_interval_list) * length(perturb_strength_list))

## lists to collect expt1rated data
experiment_params <- list()
latent_curves <- list()
for(idx in seq_len(N_CURVES)){
  
  (bw <-  runif(1, min = param_range[1, 'bandwidth'], max = param_range[2, 'bandwidth']))
  (sigma_f <- runif(1, min = param_range[1, 'sigma_f'], max = param_range[2, 'sigma_f']))
  
  t0list <- integer(length = N_EXPERIMENTS)
  delta_strength <- numeric(length = N_EXPERIMENTS)
  delta_extent <- numeric(length = N_EXPERIMENTS)
  perturbed_list <- list() # list of perturbed curves for this parameter combination
  
  ## evaluate squared-exponential kernel
  cov_kernel <- sigma_f^2 * exp(-0.5 * r_squared / (bw^2) ) + 1e-6 * diag(T) # add nugget term to make positive definite
  u_baseline <- mvrnorm(n=1, mu = rep(0, T), Sigma = cov_kernel)
  
  counter = 1 # set counter variable for experiments
  for(multiplier in perturb_strength_list){
    # define perturbation strength
    w_max <- multiplier * sigma_f # max strength of deviation is a multiple of signal noise
    # loop over intervals
    for(delta in perturb_interval_list){
      d = delta * T # width of interval
      t0 <- sample(T-(d+1), 1) # start point is in 1:{T-d}
      pert_idx <- seq(t0, t0 + d) # perturbation interval
      # compute noise function
      w <- numeric(length = T)
      w[pert_idx] <-  w_max * 0.5 * (1 - cos(2*pi*(pert_idx - t0) / d)) # noise function has same form in all cases
      perturbed_list[[counter]] = u_baseline + w
      
      ## save experiment variables
      t0list[counter] <- t0 # append to list
      delta_strength[counter] <- multiplier
      delta_extent[counter] <- delta
      
      counter <- counter + 1
    }
  }
  
  experiment_params[[idx]] <- list('sigma_f' = sigma_f,
                                   'bw' = bw, 
                                   'perturb_start' = t0list,
                                   'perturb_strength' = delta_strength,
                                   'perturb_extent' = delta_extent)
  
  latent_curves[[idx]] <- list('baseline' = u_baseline, 
                               'perturbed' = perturbed_list)
  
}


#--------------------------------------------------------------------
## use latent processes to simulate observed data
#--------------------------------------------------------------------
## for every curve, I store the data as an (n_sample x n_condition) x observed-time-points matrix 
## define experiment structure
n_samples <- 3 # 3 observations per condition 
condition_vector <- factor(rep(c('A', 'B'), each = n_samples)); K <- length(levels(condition_vector))
W <- model.matrix(~0 + condition_vector)
B <- tcrossprod(W)
noise_diag <- as.numeric(condition_vector)


## inference is made on downsampled data
## we then fit counts to original time points and compare with ground truth
obs_bins <- seq(1, T, by = 4) # observed bins
x_observed <- domain_x[obs_bins]; N <- length(x_observed)
x_test <- domain_x; M <- length(x_test)

## set noise levels relative to signal noise
noise_min <- param_range[1, 'sigma_f'] / 2
noise_max <- param_range[1, 'sigma_f'] / 1

data <- list()
for(CURVEIDX in seq_len(N_CURVES)){
  
  u_1 <- latent_curves[[CURVEIDX]][['baseline']][obs_bins]
  data[[CURVEIDX]] <- list()
  
  for(jj in seq_len(N_EXPERIMENTS)){
    
    u_2 <- latent_curves[[CURVEIDX]][['perturbed']][[jj]][obs_bins]
    ## draw observation variance for the two conditions
    
    sigma <- runif(2, min = noise_min, max = noise_max)
    
    # expt1rate data 
    # it cycles through elements by column
    y_1 <- matrix(u_1 + rnorm(n_samples*length(u_1), mean = 0, sd = sigma[1]), nrow = n_samples, byrow = T)
    y_2 <- matrix(u_2 + rnorm(n_samples*length(u_2), mean = 0, sd = sigma[2]), nrow = n_samples, byrow = T)
    
    ## append data to list in long format so I can append them all to make data matrix
    data[[CURVEIDX]][[jj]] <- c(t(rbind(y_1, y_2)))
  }
  
}

save(experiment_params, latent_curves, data, file = 'output/simulations.rdata')  

