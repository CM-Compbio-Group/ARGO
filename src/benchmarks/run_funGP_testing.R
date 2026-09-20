#' run DGE testing using funGP on simulated counts data 
## use single variance term instead of learning condition-specific noise
## fit to fewer test points, don't need to use all 2100 pseudotime values

## fit parameters 
source('src/R/simple-likelihood.R')

## Define the kernel to fit
KERNEL <- 'RBF'
FIT_CONDITIONAL_NOISE <- TRUE
EXPT <- 'conditional_noise' # expt variant i am trying out

## load input data for running testing
load('output/simulations/fungp_inputs.rdata')

sample_meta <- bin_meta %>% select(condition, sample) %>% distinct()

## define cov kernel exponents
## pairwise distances and dist-squared
r_x <- as.matrix(dist(x_observed)); dsqx <-  r_x^2
r_t <- as.matrix(dist(x_test)); dsqt <-  r_t^2
r_tx <- rdist::cdist(x_test, x_observed) ; dsqtx <- r_tx^2

## update noise_diag to have the same term for all conditions
ifelse(FIT_CONDITIONAL_NOISE, 
       noise_diag <- as.numeric(factor(sample_meta$condition)), # condition specific noise
       noise_diag <- rep(1, length(sample_meta$sample)) # all samples have same noise term
       ) 

n_noise <- length(unique(noise_diag)) # no. of noise parameters

## define the coregionalization matrix
W <- model.matrix(~0 + factor(sample_meta$condition))
B <- tcrossprod(W)

## call the parameter estimation script
## run gp parameter inference
## pairwise distances and dist-squared
min_bw <- min(dist(x_observed)) * 0.5
min_log_bw <- log(min_bw)
max_bw <- 3*diff(range(x_observed))
max_log_bw <- log(max_bw)


kernelParams <- list()
for(gene in genes){
  message(gene)
  y <- bin_meta %>% pull(gene)
  log_bw <- runif(1, log(.1), log(.5))
  
  # compute empirical variation in expt1 by condition group
  emp_sd <- bin_meta %>% 
              {if(FIT_CONDITIONAL_NOISE) group_by(., condition) else .} %>%
               summarize(sd_val = sd(!!sym(gene))) %>% 
               pull(sd_val)
  
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
                             KERNEL = KERNEL, # fit specific covariance-kernel to data
                             method = 'L-BFGS-B', 
                             lower = c(min_log_bw, rep(-3, 1 + n_noise)),
                             upper = c(max_log_bw,  rep(3, 1 + n_noise)),
                             control = list(maxit = 1000, trace = 2))
  
  kernelParams[[gene]] <- c(bandwidth = exp(optimResult$par[1]), 
                            sigma_f = exp(optimResult$par[2]), 
                            sigma_n = exp(optimResult$par[3:(3 + n_noise - 1)]))
  
  
}

params <- do.call('rbind', kernelParams)
colnames(params) <- c('bandwidth', 'sigma_f', paste('sigma_n', 1:n_noise, sep = ''))
pheatmap(params)

## save inferred params to file
write_rds(params,glue('output/simulations/funGP_params_{EXPT}.rds'))

##-------------------------------------------------------------------
## visualize posterior fits for some experiments
##-------------------------------------------------------------------
source('src/R/simple_posterior_fits.R')
genes_to_plot <- sample(effect_genes, 8)
plots <- list()
for(gene in genes_to_plot){
  message(gene)
  ## pass corresponding gp parameters
  plots[[gene]] <- simple_posterior_fit(gene)
}

ggsave(wrap_plots(plots, ncol = 4)+ 
         plot_layout(guides = "collect"), 
       filename = file.path('figs', glue('funGP_sample_fits.pdf')),
       width = 15, height = 8)

#===================================================================
## run differential testing with funGP
#===================================================================
source('src/R/simple_contrast_test_with_funGP.R')
source('src/R/generate_cc_matrix.R')
source('src/R/compute_analytical_pvalue.R')
source('src/R/derive_KLT_CI_bands.R')

## control variables
KLT_CONF_BAND <- TRUE
CONTRAST_VECTOR <- c(-1, 1)

funGP_fits <- vector("list", G)
for(gene in genes){
  message(gene)
  funGP_fits[[gene]]  <- simple_contrast_test_with_funGP(contrast_vector = CONTRAST_VECTOR, gene = gene)
}

## save diff testing results to file
write_rds(funGP_fits, glue('output/simulations/funGP_res_{EXPT}.rds'))

#===================================================================
## gather results and visualize confidence bands
#===================================================================
source('src/R/plot_conf_bands.R')
source('src/R/compute_exceedence_metrics.R')
source('src/R/compute_exceedence_statistics.R')

stat_binary <- compute_stat_binary(funGP_fits, M)
binary_colors <- c("1" = "red3","0" = "lightgrey", "-1" = "navy")
pheatmap(stat_binary, cluster_cols = F)

pdf(glue('figs/funGP_binary_heatmap_{EXPT}.pdf'), width = 8)
plot(Heatmap(stat_binary, 
             cluster_columns = F, cluster_rows = T,
             name = 'binary \nstatistic',
             show_row_names = T,
             row_names_gp = gpar(fontsize=4),
             column_names_gp = gpar(fontsize=4),
             col = binary_colors,
             heatmap_legend_param = list(color_bar = "discrete"), use_raster = F
))
dev.off()

#-----------------------------------------------------
## conf bands for selected experiments
#-----------------------------------------------------
conf_band_plots <- lapply(genes_to_plot, function(gene) {
  plot_CI_bands(gene, x = x_test, 
                band = funGP_fits[[gene]][['band']],
                stat = funGP_fits[[gene]][['statistic']], 
                alpha = 95)
  }
)

ggsave(wrap_plots(conf_band_plots, ncol = 4) + 
         plot_layout(guides = "collect"), 
       filename = file.path('figs', glue('funGP_conf_bands_{EXPT}.pdf')),
       width=18,height=8)

#-----------------------------------------------------
## summarize experiment results
#-----------------------------------------------------
## compute exceedence metrics
# exc.metrics <- do.call('rbind', lapply(funGP_fits, compute_exceedence_metrics, x_test))

exc.metrics <- list()
for(gene in genes){
  exc.metrics[[gene]] <- compute_exceedence_metrics(funGP_fits[[gene]], x_test)
}
exc.metrics <- do.call('rbind', exc.metrics)

## aggregate results from functional testing
## append gene conditions
funGP.df <- funGP_fits %>%
  map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'gene') %>%
  bind_cols(exc.metrics) %>% 
  mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH'),
         net_auec = auec_up - auec_down,
         net_poe = poe_up - poe_down) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())

funGP.df

p <- get_ranked_dge_plot(funGP.df, method = 'funGP', markers = effect_genes)
print(p)

funGP.df %>% filter(fdr < 0.05) %>% with(table(deviation > 0))
funGP.df %>% filter(fdr < 0.05 & max_deviation > 0) %>% with(table(gene %in% null_genes)) ## how many control genes are being called differential

boxplot(funGP.df$deviation~ as.numeric(funGP.df$gene %in% null_genes))
write_rds(funGP.df, glue('output/simulations/funGP_res_{EXPT}.rds'))

#===================================================================
## compare recovered effects with ground truth
#===================================================================
## compare performance on effect genes
res_effect <- filter(funGP.df, gene %in% effect_genes) %>% 
  arrange(match(gene, effect_info$gene)) %>% right_join(effect_info)

sum(res_effect$fdr < .05)
filter(res_effect, fdr > .05)
plot(res_effect$index~ res_effect$center)
plot(res_effect$deviation~res_effect$amplitude)
plot(res_effect$deviation~res_effect$bandwidth)
table(res_effect$deviation > 0, res_effect$sign)



