#' read scDesign parameter estimates for dataset and perform perturbation experiments
source('src/benchmarks/00_load_packages_and_functions.R')

start_time <- Sys.time()

dataset <- 'PANCREAS'
load(glue('output/simulations/{dataset}_scdesign.rdata'))

t_ <- PANCREAS_data$dat$pseudotime # pseudotime values for all cells
ord <- order(t_)

N_EXPERIMENTS <- 10
N_SAMPLES <- 3 # no. of observations per condition
EXPT_params <- vector("list", length = N_EXPERIMENTS)
RES_list <- vector("list", length = N_EXPERIMENTS)
fdr_thresh <- 0.05

## iterte through a number of experiments
for(iter in seq_len(N_EXPERIMENTS)){
  set.seed(148 + iter)
  
  ## sample genes to include in experiment
  genes <- sample(colnames(PANCREAS_data$count_mat), size = 500, replace = F); 
  G <- length(genes)
  ## extract the mean parameters for these genes
  mu_B <- PANCREAS_para$mean_mat[, genes]
  
  #====================================
  ## decide which genes to perturb
  #====================================
  effect_genes <- sample(genes, size = G/5, replace = F)
  null_genes <- setdiff(genes, effect_genes) # useful to know false positives in differential testing
  
  lfc_mat <- matrix(
    0,
    nrow = nrow(mu_B),
    ncol = G,
    dimnames = dimnames(mu_B)
  )
  
  #====================================
  # compute effects
  #====================================
  effect_info <- data.frame(
    gene = effect_genes,
    amplitude = sample(c(0.1, .25, .5, .75, 1, 1.5), size = length(effect_genes), replace = T), # sample discretely
    sign = sample(c(-1, 1), size = length(effect_genes), replace = T),
    center = runif(length(effect_genes), 0.2, 0.8),
    bandwidth = sample(
      c(0.025, 0.05, 0.075, 0.10, 0.15),
      length(effect_genes),
      replace = TRUE
    )
  )
  
  ## inject perturbations on mean parameter
  for (j in seq_len(nrow(effect_info))) {
    g <- effect_info$gene[j]
    
    lfc_mat[, g] <-
      effect_info$sign[j] * effect_info$amplitude[j] *
      exp(
        -(t_ - effect_info$center[j])^2 /
          (2 * effect_info$bandwidth[j]^2)
      )
  }
  
  #====================================
  ## visualize ground truth effects
  #====================================
  # pdf(glue('figs/{dataset}_{iter}_mean_perturbations.pdf'), width = 12)
  # pheatmap(lfc_mat[ord, effect_genes], cluster_rows = F, 
  #          show_rownames = F, fontsize_col = 9, breaks = c(-1,0, 1), 
  #          use_raster=F, name = 'LFC', )
  # dev.off()
  
  #===================================================================
  ## simulate counts data
  #===================================================================
  ## compute the new mean parameter
  mu_B <- mu_B * exp(lfc_mat)
  copula_list <- PANCREAS_copula$copula_list
  copula_list$`1` <- copula_list$`1`[genes, genes]
  
  generate_counts <- function(mu){
    simu_new(sce = pancreas_sce[genes, ],
             mean_mat = mu,
             sigma_mat = PANCREAS_para$sigma_mat[, genes],
             zero_mat = PANCREAS_para$zero_mat[, genes],
             quantile_mat = NULL,
             copula_list = copula_list,
             n_cores = 20,
             family_use = "nb",
             input_data = PANCREAS_data$dat,
             new_covariate = PANCREAS_data$newCovariate,
             filtered_gene = PANCREAS_data$filtered_gene)
  }
  
  counts_A <- vector("list", N_SAMPLES)
  counts_B <- vector("list", N_SAMPLES)
  
  for(r in seq_len(N_SAMPLES)){
    set.seed(135 + iter + r)
    counts_A[[r]] <- generate_counts(PANCREAS_para$mean_mat[, genes])
    counts_B[[r]] <- generate_counts(mu_B)
  }
  
  
  ## save run parameters to file
  EXPT_params[[iter]] <- list(
    run_iter = iter,
    sampled_genes = genes, 
    effect_info = effect_info,
    effect_genes = effect_genes
  )
  
  
  sim_counts <- do.call(cbind, c(counts_A, counts_B))
  ## append unique sample name for each datapoint
  sim_meta <- data.frame(condition = factor(rep(c('condition_A', 'condition_B'), 
                                                each = N_SAMPLES * ncol(pancreas_sce))),
                         replicate = rep(factor(rep(paste('replicate', 1:N_SAMPLES, sep = '_'), 
                                                    each = ncol(pancreas_sce))), 
                                         2), 
                         sample = rep(factor(paste('sample', LETTERS[1:6], sep = '-')), 
                                      each = ncol(pancreas_sce)),
                         pseudotime = rep(PANCREAS_data$dat$pseudotime, 2*N_SAMPLES)) %>% as_tibble()
  
  
  #===================================================================
  ## prepare input data for differential testing
  #===================================================================
  sim_meta <-  sim_meta %>% 
    mutate(bin_observed = factor(ntile(pseudotime, n = 25)), ## 25 points where data is observed
           bin_test = factor(ntile(pseudotime, n = 100)), # 100 points to test fits at
           sample_bin = paste(sample, bin_observed, sep = '_'),
           condition_bin = paste(condition, bin_observed, sep = '_')) %>% 
            mutate(condition_bin = factor(condition_bin, 
                                  levels =  gtools::mixedsort(unique(condition_bin))),
                  sample_bin = factor(sample_bin, 
                                  levels = gtools::mixedsort(unique(sample_bin)))) 
  
  ## append coordinates for each bin
  sim_meta <- right_join(sim_meta, sim_meta %>% group_by(bin_observed) %>% 
                           summarise(coord_x = mean(pseudotime)), by = 'bin_observed')
  sim_meta <- right_join(sim_meta, sim_meta %>% group_by(bin_test) %>% 
                           summarise(test_x = mean(pseudotime)), by = 'bin_test')
  sim_meta
  
  ## prepare bin level metadata
  bin_level_meta <- sim_meta %>% dplyr::count(condition, replicate, sample, 
                                        condition_bin, sample_bin,
                                        coord_x, bin_observed,
                                        name = 'n_cells') %>% 
    distinct() %>% 
    arrange(sample_bin)
  
  ## sample level information
  sample_meta <- bin_level_meta %>% select(condition, sample) %>% distinct()
  
  ## convert counts to pseudobulk
  binned_counts <- rowsum(t(sim_counts), group = sim_meta$sample_bin)
  all(rownames(binned_counts) == bin_level_meta$sample_bin)
  
  #===================================================================
  ## prepare deseq object with pooled bin-level counts data
  #===================================================================
  dds <- DESeqDataSetFromMatrix(countData = t(binned_counts),
                                colData = bin_level_meta,
                                design = ~ 0 + condition_bin)
  
  dds <- estimateSizeFactors(dds)
  ## perform counts normalization
  vsd <- varianceStabilizingTransformation(dds, blind = F)
  
  #===================================================================
  ## append the counts to bin meta to prepare input for funGP
  #===================================================================
  ## append the counts data to the bin_meta object
  all(colnames(vsd) == bin_level_meta$sample_bin)
  ## append centered log-normalized counts data
  bin_meta <- bin_level_meta %>%
    # mutate(bin = as.character(bin)) %>%
    bind_cols(scale(t(assay(vsd)), center = T, scale = F))
  
  ## get variables for fitting continuous-time models
  x_observed <- sort(unique(bin_meta$coord_x)); N <- length(x_observed)
  x_test <- sort(unique(sim_meta$test_x)); M <- length(x_test)
  
  #===================================================================
  ## call differential testing methods
  #===================================================================
  message('running diff testing with baseline methods: \n')
  dge.res.list <- run_baseline_methods(dds = dds, binned_counts =  binned_counts, 
                                   bin_level_meta = bin_level_meta, sample_meta = sample_meta, 
                                   N = N)
  ## call source codes to fit argo and run testing
  message('running diff testing with argo')
  dge.res.list[['argo_shared_noise']] <- run_argo(FIT_CONDITIONAL_NOISE = 'FALSE')
  dge.res.list[['argo_conditional_noise']] <- run_argo(FIT_CONDITIONAL_NOISE = 'TRUE')
  
  
  sapply(dge.res.list, function(df) sum(df[['fdr']] < fdr_thresh))
  eval.list <- lapply(dge.res.list, compute_eval_metrics, effect_genes, null_genes)
  
  
  ## append evaluation metrics to EXPT list
  RES_list[[iter]] <- list(dge.res = dge.res.list, 
                          eval.metrics = do.call('rbind', lapply(eval.list, '[[', 'metrics'))
                      )
  
  
  SX <- glue('pancreatic_endocrinogenesis-{iter}')
  cols <- c( 'burlywood2', 'lightskyblue2', 'lightgreen', 'lightpink', 'lightsalmon', 'plum3', 'orchid3')
  
  pr.curves <- lapply(eval.list, '[[', 'probj')
  
  pdf(glue('figs/AUPRC_{SX}.pdf'), width = 8, height = 8)
  par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
  plot(pr.curves[['deseq']], col = adjustcolor(cols[1], alpha.f = 1), 
       main = glue("{SX} \n Precision-Recall Curve"), 
       auc.main = F, lty=2, lwd = 2)
  ## plot other curves
  for(ii in seq(2,7)){
    method = names(pr.curves)[ii]
    plot(pr.curves[[method]], col = adjustcolor(cols[ii + 1], alpha.f = .5), lty = 2, lwd = 2, add = TRUE)
  }
  
  ## make legends
  lgx <- glue("DESeq2 (AUPRC = ", round(pr.curves[['deseq']]$auc.integral, 3), ")")
  for(method in names(pr.curves)[2:7]){
    lgx <- c(lgx, 
             glue("{method} (AUPRC = {round(pr.curves[[method]]$auc.integral, 3)}) ")
    )
  }
  legend("bottomleft", cex=.8,
         legend = lgx,
         col = cols, 
         lty = 2,
         lwd = 1,             # Line width
         bty = "n")           # Removes the box border around the legend (optional)
  
  dev.off()
  
  
}

end_time <- Sys.time()
print(end_time - start_time)

save(RES_list, EXPT_params, file = glue('output/simulations/{dataset}_results.rdata'))
