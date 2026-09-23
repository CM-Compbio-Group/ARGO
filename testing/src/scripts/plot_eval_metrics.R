## gather all the results and make boxplots comparing performance of the different methods

load(glue('output/simulations/{dataset}_results.rdata'))

df <- EXPT_params[[1]]$effect_info
table(df$amplitude, df$bandwidth)

#===================================================================
## compute evaluation metrics for different experimental parameters
#===================================================================
## each iteration has a different set of true and null genes
## each effect has a different set of ground-truth genes

## gather effect info across all experiments and jointly filter

for(iter in seq_len(N_EXPERIMENTS)){
  message(iter)
  dge.list <- RES_list[[iter]]$dge.res
  null_genes <- setdiff(EXPT_params[[iter]]$sampled_genes, EXPT_params[[iter]]$effect_genes) # this list doesn't change by effect
  
}

#===================================================================
## plot overall performance metrics in single figure
#===================================================================
eval.df <- do.call('rbind', lapply(RES_list, '[[', 'eval.metrics')) 
method <- rownames(eval.df)

eval.df <- as_tibble(eval.df) %>% 
  mutate('method' = method,
    'iter' = rep(seq(10), each = 7))
  

library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Reshape the data from wide to long format
eval_long <- eval.df %>%
  mutate(method = as.factor(method)) %>%
  pivot_longer(
    cols = c(sensitivity, specificity, precision, empirical_FDR, accuracy, F1, AUPRC),
    names_to = "metric",
    values_to = "value"
  )

# 2. Generate the plot grid
p <- ggplot(eval_long, aes(x = method, y = value, color = method)) +
  # Scales the 10 iterations tightly within each of the 7 method columns
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  # Adds a clear center-point for the mean
  stat_summary(fun = mean, geom = "point", color = "black", size = 2, shape = 18) +
  # Splitting into metric subpanels
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") + 
  labs(
    x = "Method",
    y = "Score",
    title = "Metric Performance Across 10 Iterations"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(vjust = 1, hjust = 1, size = 12),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing.x = unit(1.5, "lines"),  # Horizontal space between columns
    panel.spacing.y = unit(2.5, "lines")   # Increased vertical space between rows
  )

ggsave(plot = p, filename = glue('figs/eval_metrics_{dataset}.pdf'), width = 12, height = 6,device = 'pdf')
