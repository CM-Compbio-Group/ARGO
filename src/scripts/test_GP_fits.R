# source('src/R/posterior_predictive_plots.R')
# source('src/R/nll_and_gradients.R')
# source('src/R/plot_conf_bands.R')

message(gene)
## abundance variable to be modeled
y <- bin_meta %>% pull(gene)
## draw bandwidth parameter
log_bw <- runif(1, log(.2), log(.5))
## sigma_f is set to mean SD of y by condition
emp_sd <- bin_meta %>% group_by(age) %>%
  summarize(sd_val = sd(!!sym(gene))) %>%
  pull(sd_val)
log_sigma_f <-  log(mean(emp_sd) * sqrt(0.7)) ## 70% variance to signal
# log_sigma_n <- rep(log(.15* sqrt(.3)), 3)
log_sigma_n <- log(emp_sd * sqrt(0.3)) ## 30% variance to noise
parInt <- c(log_bw, log_sigma_f, log_sigma_n)

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
                           d_squared = dsqx,
                           noise_diag = noise_diag,
                           slides_per_group = slides_per_age,
                           B = B, 
                           lambda = 1,
                           method = 'L-BFGS-B', 
                           lower = c(min_log_bw, -3, rep(-5, 3)),
                           upper = c(max_log_bw, 3, rep(5, 3)),
                           control = list(maxit = 1000, trace = 2))

estimatedParams = list(bandwidth = exp(optimResult$par[1]), 
                       sigma_f = exp(optimResult$par[2]), 
                       sigma_n = exp(optimResult$par[3:5]))

## extract and store the parameters 
print(par <- unlist(estimatedParams))

posterior_predictive_plots(cell = gene, y = y,
                          par = par, 
                          x_observed = x_observed,
                          x_test = x_test,
                          P_obs = P_obs,
                          dsqx = dsqx,
                          dsqt = dsqt,
                          B = B, 
                          N = N, 
                          M = M,
                          slides_per_group = slides_per_age,
                          noise_diag = noise_diag,
                          bin_meta = bin_meta)

## get test results
res <- contrast_test_with_funGP(contrast_vector = c(-1, 0, 1), 
                         y = y, 
                         P_obs = P_obs,
                         par = par, 
                         x_observed = x_observed,
                         x_test = x_test,
                         N = N,
                         M = M,
                         dsqx = dsqx,
                         dsqt = dsqt,
                         W = W,
                         B = B, 
                         noise_diag = noise_diag, 
                         return_conf_band = T,
                         alpha = 0.95
)


plot_CI_bands(gene, x = x_test, 
              band = res[['band']],
              stat = res[['statistic']], 
              alpha = 95)
