#' read the DGE testing evaluation metrics and make boxplots of AUPRC values and F1 scores for each method
#' and save PDF files to folder 'testing/figs/'

suppressMessages({
  source('testing/src/scripts/00_load_packages_and_functions.R')
})

expt_df <- vroom::vroom('testing/data/temp/experiment_params.csv')
eval.metrics <- read_rds('testing/output/diff_testing_eval_metrics.rds')
metrics <- do.call(rbind, eval.metrics)
method <- rownames(metrics)
eval_long <- metrics %>% as_tibble() %>% 
  mutate(method = method, 
         Experiment = rep(seq(nrow(expt_df)), each = 7)) %>%
  mutate(method = factor(method, levels = c('DESeq2', 'edgeR', 'Limma', 'tradeSeq', 'Lamian', 'argo_shared_noise', 'argo_conditional_noise'))) %>%
  pivot_longer(
    cols = c(sensitivity, specificity, precision, empirical_FDR, accuracy, F1, AUPRC),
    names_to = "metric",
    values_to = "value"
  )

eval_long <- merge(eval_long, expt_df, by = 'Experiment') %>% as_tibble()
eval_long

for(metrx in c('AUPRC', 'F1')){
  message(metrx)
  p <- ggplot(eval_long %>% filter(metric == metrx), aes(x = method, y = value, color = method)) +
    # Scales the 10 iterations tightly within each of the 7 method columns
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) + 
    # Adds a clear center-point for the mean
    stat_summary(fun = mean, geom = "point", color = "black", size=3, pch=18) +
    # Splitting into metric subpanels
    facet_wrap(~ amplitude * bandwidth, scales = "free_y", ncol = 5, 
               labeller = \(x) label_both(x, multi_line = TRUE, sep = ": ")) +
    theme_minimal(base_size = 12) +
    scale_color_brewer(palette = "Set2") + 
    labs(
      x = "",
      y = metrx,
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8),
      axis.text.y = element_text(vjust = 1, hjust = 1, size = 8),
      legend.position = "none",
      strip.text = element_text(face = "plain", size = 6,
                                margin = margin(t = .5, b = .5, r = .5, l = .5)
                                ),
      panel.spacing.x = unit(0.1, "lines"),  # Horizontal space between columns
      panel.spacing.y = unit(0.1, "lines"),   # Increased vertical space between rows
      plot.margin = margin(t = 1, r = 2, b = 1, l = 2, unit = "pt") 
    ) + 
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
  
  
  ggsave(plot = p, filename = glue('testing/figs/{metrx}_evaluation.pdf'), 
         device = 'pdf', width = 14, height = 11.69, units='in')
  
}
