#' **02_visualize_fitted_GPs** : use the fitted GP parameters to visualize functional curve fit for each gene along isodepth
#' @params kernel_params : fitted GP parameters per gene
#' @params genes_to_plot : genes of interest to visualize posterir fits for
#' returns list of fitted plots for the genes


## make the posterior predictive plots for given list of genes
source('src/R/posterior_predictive_plots.R')
visualize_fitted_GP <- function(genes_to_plot, params, FILTER_TO_CONTRAST = F){
  plots <- list()
  for(gene in genes_to_plot){
    message(gene)
    y <- bin_meta %>% pull(gene)
    ## pass corresponding gp parameters
    plots[[gene]] <- posterior_predictive_plots(gene = gene, y = y,
                                                par = params[gene, ], 
                                                x_observed = x_observed,
                                                x_test = x_test,
                                                P_obs = P_obs,
                                                r_x = r_x,
                                                r_t = r_t,
                                                r_tx = r_tx,
                                                dsqx = dsqx,
                                                dsqt = dsqt,
                                                dsqtx = dsqtx,
                                                B = B, 
                                                N = N, 
                                                M = M,
                                                slides_per_group = slides_per_age,
                                                noise_diag = noise_diag,
                                                KERNEL = KERNEL, 
                                                fit_donor = FIT_DONOR_SWITCH,
                                                cov_donor = cov_donor,
                                                meta = bin_meta,
                                                slide_df = slide_df,
                                                FILTER = FILTER_TO_CONTRAST)
  }
  return(plots)
}



# plots <- visualize_fitted_GP(sample(genes, 5), params)
# 
# ggsave(wrap_plots(plots, ncol = 5) + 
#              plot_layout(guides = "collect"), 
#            filename = file.path('figs', glue('{KERNEL}_fit_sample.pdf')),
#            width=18,height=4)
