## fit parameters 

## define cov kernel exponents
## pairwise distances and dist-squared
r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2

## call the parameter estimation script
## run gp parameter inference
## pairwise distances and dist-squared
min_bw <- min(dist(x_observed)) * 0.5
min_log_bw <- log(min_bw)
max_bw <- 3*diff(range(x_observed))
max_log_bw <- log(max_bw)


y <- simx.data %>% pull(gene)
## draw bandwidth parameter
log_bw <- runif(1, log(.2), log(.5))
# compute empirical variation in expt1 by condition group
emp_sd <- simx.data %>% group_by(condition) %>% summarize(sd_val = sd(!!sym(gene))) %>% pull(sd_val)
log_sigma_f <-  log(mean(emp_sd) * sqrt(0.7)) ## 70% variance to signal
log_sigma_n <- log(emp_sd * sqrt(0.3)) ## 30% variance to noise
parInt <- c(log_bw, log_sigma_f, log_sigma_n)

optimResult = stats::optim(par = parInt, 
                           fn = computeNegLogLikelihood, 
                           gr = computeNegLogLikelihoodGradient,
                           x = x_observed,
                           y = y,
                           N = N, 
                           r = r_x,
                           r_squared = dsqx,
                           noise_diag = noise_diag,
                           B = B,
                           lambda = 1, # regularization for noise nuggets
                           kernel = KERNEL, # fit specific covariance-kernel to data
                           method = 'L-BFGS-B', 
                           lower = c(min_log_bw, rep(-3, 3)),
                           upper = c(max_log_bw,  rep(3, 3)),
                           control = list(maxit = 1000, trace = 2))

params = c(bandwidth = exp(optimResult$par[1]), 
           sigma_f = exp(optimResult$par[2]), 
           sigma_n = exp(optimResult$par[3:4]))

params

## compare fitted params with original values
experiment_params[[1]][1:2] ## fitted bandwidth is much lower



