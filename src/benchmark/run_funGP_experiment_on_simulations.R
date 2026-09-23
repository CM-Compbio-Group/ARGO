#' we now have N_CURVES x n_perturbations experiments
#' for each experiment, we can draw samples y_i = u_i + \epsilon_i with condition-specific observation noise
source('src/R/simple-likelihood.R')

load('output/simulations.rdata')

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

funGP.res <- list()
## perform inference on one set of experiments
for(CURVEIDX in seq_len(N_CURVES)){
  
  message(CURVEIDX)
  ## visualize set of genrated perturbation curves
  U_perturb <- do.call('rbind', latent_curves[[CURVEIDX]]$perturbed)
  
  pdf(glue('figs/perturbed_processes_{CURVEIDX}.pdf'), width = 8)
  matplot(domain_x, t(U_perturb), type = "b", lty = 1, pch = 19, cex=.25, 
          lwd = .5, col = 1:nrow(U_perturb), ylim = range(c(U_perturb)),
          xlab = "Columns", ylab = "Values", main = "perturbed processes")
  matplot(domain_x, latent_curves[[CURVEIDX]]$baseline, add = T, pch = 20, cex=.5, lwd=1, col = 'red4', lty = 2)
  dev.off()
  
  
  ## create the data matrix 
  simx.data <- tibble(coord_x = rep(x_observed, K*n_samples),
                      SampleID = rep(rep(paste('Sample', 1:n_samples, sep='-'), each = N),
                                     K),
                      condition = rep(levels(condition_vector), each = n_samples * N))
  
  ## append simulated counts data
  simx.data <- bind_cols(simx.data,  do.call('cbind', data[[CURVEIDX]]))
  colnames(simx.data)[4:ncol(simx.data)] <- paste('Experiment', 1:N_EXPERIMENTS, sep = '_')
  
  # ggplot(simx.data, aes(x = coord_x, y = `Experiment_18`, col = condition, shape = SampleID)) + 
  #   geom_point(alpha = 0.5, cex = 3.5)
  
  
  ## independently infer parameter values from each set of observations
  kernelParams <- list()
  for(idx in seq(N_EXPERIMENTS)){
    expt <- paste('Experiment', idx, sep = '_')
    y <- simx.data %>% pull(expt)
    ## draw bandwidth parameter
    log_bw <- runif(1, log(.2), log(.5))
    # compute empirical variation in expt1 by condition group
    emp_sd <- simx.data %>% group_by(condition) %>% summarize(sd_val = sd(!!sym(expt))) %>% pull(sd_val)
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
    
    kernelParams[[expt]] <- c(bandwidth = exp(optimResult$par[1]), 
                             sigma_f = exp(optimResult$par[2]), 
                             sigma_n = exp(optimResult$par[3:4]))
    
  }
  
  params <- do.call('rbind', kernelParams)
  colnames(params) <- c('bandwidth', 'sigma_f', 'sigma_x1', 'sigma_x2')
  
  ## save inferred params to file
  write_rds(params,glue('output/Expt_{CURVEIDX}_params.rds'))
  ## compare with experiment settings
  # unlist(experiment_params[[CURVEIDX]][1:2])
  
  ##-------------------------------------------------------------------
  ## visualize posterior fits for some experiments
  ##-------------------------------------------------------------------
  source('src/R/simple_posterior_fits.R')
  
  expt_to_plot <- sort(sample(N_EXPERIMENTS, 8))
  plots <- list()
  for(idx in expt_to_plot){
    message(idx)
    expt <- paste('Experiment', idx, sep = '_')
    ## pass corresponding gp parameters
    plots[[expt]] <- simple_posterior_fit(expt)
  }
  
  ggsave(wrap_plots(plots, ncol = 4)+ 
        plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('Expt_{CURVEIDX}_sample_fits.pdf')),
        width = 15, height = 8)
  
  
  #--------------------------------------------------------------------
  ## perform inference over all experiments and fit confidence bands
  #--------------------------------------------------------------------
  source('src/R/generate_cc_matrix.R')
  source('src/R/compute_analytical_pvalue.R')
  source('src/R/derive_KLT_CI_bands.R')
  source('src/R/simple_contrast_test_with_funGP.R')
  
  ## control variables
  KLT_CONF_BAND <- TRUE
  CONTRAST_VECTOR <- c(-1, 1)
  
  funGP_fits <- list()
  for(idx in seq(N_EXPERIMENTS)){
    expt <- paste('Experiment', idx, sep = '_')
    message(expt)
    funGP_fits[[expt]]  <- simple_contrast_test_with_funGP(contrast_vector = CONTRAST_VECTOR, expt = expt)
  }
  
  ## save diff testing results to file
  write_rds(funGP_fits, glue('output/Experiment_{CURVEIDX}_res.rds'))
  
  #-----------------------------------------------------
  ## gather results and visualize conf bands
  ## visualize outputs of fungp tests
  #-----------------------------------------------------
  source('src/R/plot_conf_bands.R')
  source('src/R/compute_exceedence_metrics.R')
  source('src/R/compute_exceedence_statistics.R')
  
  stat_binary <- compute_stat_binary(funGP_fits, M)
  binary_colors <- c("1" = "red3","0" = "lightgrey", "-1" = "navy")
  pdf(glue('figs/expt_{CURVEIDX}_binary_heatmap.pdf'), width = 8)
  plot(Heatmap(stat_binary, 
          cluster_columns = F, cluster_rows = F,
          name = 'binary \nstatistic',
          show_row_names = T,
          row_names_gp = gpar(fontsize=8),
          column_names_gp = gpar(fontsize=4),
          col = binary_colors,
          heatmap_legend_param = list(color_bar = "discrete")
  ))
  dev.off()
  
  # stat_lfc <- compute_stat_continuous(funGP_fits, M)
  # Heatmap(stat_lfc, cluster_columns = F, cluster_rows = F,
  #         name = 'FC \nstatistic',
  #         row_names_gp = gpar(fontsize=8),
  #         column_names_gp = gpar(fontsize=4))
  
  #-----------------------------------------------------
  ## conf bands for selected experiments
  #-----------------------------------------------------
  conf_band_plots <- lapply(expt_to_plot, function(idx) {
      expt <- paste('Experiment', idx, sep = '_'); 
      plot_CI_bands(expt, x = x_test, 
                    band = funGP_fits[[expt]][['band']],
                    stat = funGP_fits[[expt]][['statistic']], 
                    alpha = 95)
    }
  )
  
  ggsave(wrap_plots(conf_band_plots, ncol = 4) + 
          plot_layout(guides = "collect"), 
          filename = file.path('figs', glue('Expt_{CURVEIDX}_conf_bands.pdf')),
          width=18,height=8)
  
  #-----------------------------------------------------
  ## summarize experiment results
  #-----------------------------------------------------
  ## compute exceedence metrics
  exc.metrics <- do.call('rbind', lapply(funGP_fits, compute_exceedence_metrics, x_test))
  
  ## aggregate results from functional testing
  ## append expt conditions
  funGP.df <- funGP_fits %>%
    map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'expt') %>%
    mutate(perturb_start = experiment_params[[CURVEIDX]]$perturb_start, 
           perturb_strength = experiment_params[[CURVEIDX]]$perturb_strength,
           perturb_extent = experiment_params[[CURVEIDX]]$perturb_extent) %>%
    bind_cols(exc.metrics) %>% 
    mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
    mutate(stats = neglogp * sign(deviation), 
           net_auec = auec_up - auec_down,
           net_poe = poe_up - poe_down) %>%
    arrange(desc(neglogp)) %>%
    mutate(Rank = row_number())
  
  
  
  # plot(funGP.df$neglogp~ funGP.df$perturb_strength, pch = 20, cex = 1.5, col = adjustcolor('red', alpha=.5))
  # plot(funGP.df$neglogp~ funGP.df$perturb_extent, pch = 20, cex = 1.5, col = adjustcolor('red', alpha=.5))  
  
  funGP.res[[CURVEIDX]] <- funGP.df
}

write_rds(funGP.res, file = glue('output/funGP_simulated_testint_results.rds'))