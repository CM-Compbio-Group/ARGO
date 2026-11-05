# sid <- 29
# 
# # load simulation data
# dat <- readRDS(str_interp('${datadir}/PANCREAS_expt_${sid}.rds'))
# effect_genes <- dat$effect_genes
# amplitude <- dat$effect_info$amplitude[1]
# bandwidth <- dat$effect_info$bandwidth[1]
#   
# # load SPADE prediction
# for (rep in seq(1, 4)){
#   if (file.exists(str_interp("${spade_dir}/PANCREAS_expt_${sid}_replicate${rep}.tsv"))) {
#     spade_res <- read.table(str_interp("${spade_dir}/PANCREAS_expt_${sid}_replicate${rep}.tsv"))
#     pos_scores <- -log10( spade_res[spade_res$geneid %in% effect_genes, 'Adjust.Pvalue'] + 1e-50 )
#     neg_scores <- -log10( spade_res[!spade_res$geneid %in% effect_genes, 'Adjust.Pvalue'] + 1e-50 )
#     pr_curve_result <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = TRUE)
#     # store AUPRC value for this sid and rep
#     auprc <- pr_curve_result$auc.integral
#     spade_auprc <- rbind(
#       spade_auprc,
#       data.frame(sid = sid, rep = rep, amplitude = amplitude, bandwidth = bandwidth, auprc = auprc)
#     )
#   }
# }


spade_auprc %>%
  group_by(amplitude, bandwidth) %>%
  summarise(across(auprc, \(x) mean(x, na.rm = TRUE)), .groups = 'drop')

spade <- spade_auprc %>% group_by(sid, amplitude, bandwidth) %>% 
  summarise(across(auprc, \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(method = 'SPADE') %>%  dplyr::rename(Experiment = sid)

despace_auprc %>% mutate(method = 'DESpace')

despace_auprc %>%
  group_by(amplitude, bandwidth) %>%
  summarise(across(auprc, \(x) mean(x, na.rm = TRUE)), .groups = 'drop')


eval.metrics <- read_rds(glue('testing/output/diff_testing_N_25_p_500_eval_metrics.rds'))
metrics <- do.call(rbind, eval.metrics)

metrics <- metrics %>% as_tibble(rownames = "method") %>%
  select(method, auprc = AUPRC) %>%
  mutate( method = factor(method, levels = unique(method)), 
        Experiment = rep(seq(nrow(expt_df)), each = 6)) %>% 
        mutate(method = fct_recode(method, 
                             'ARGO' = 'argo_conditional_noise')) %>% 
  right_join(expt_df %>% select(-iter), by = 'Experiment')


metrics <- metrics %>% bind_rows(despace_auprc %>% 
                     mutate(method = 'DESpace') %>% 
                     dplyr::rename(Experiment = sid)
                     ) %>% 
                    bind_rows(spade) %>%
                    arrange(Experiment)

metrics <- metrics %>% right_join(rel_length, by = 'bandwidth') %>%
              mutate(across(`relative length`, ~ round(., 3))) %>% 
              dplyr::rename(`signal strength` = amplitude) %>% 
  dplyr::rename(RL = `relative length`,
                SS = `signal strength`)  %>% 
  mutate(RL = factor(RL), 
         SS = factor(SS)) 



method_cols <- c(
  "ARGO"     = "orchid4",
  "DESeq2"   = "#bebada",
  "edgeR"    = "#8dd3c7",
  "Limma"    = "lightskyblue2",
  "Lamian"   = "lightgreen",
  "tradeSeq" = "lightpink2", 
  'DESpace' = 'gold2', 
  'SPADE' = 'maroon3'
)

metrics <- metrics %>% mutate(method = fct_relevel(method, 
                                       "SPADE", "DESeq2", "edgeR", "Limma", 
                                       "DESpace", "tradeSeq", "Lamian", "ARGO"))


write_rds(metrics, 'testing/output/auprc_metrics_w_spade.rds')

## drop spade, rename despace 
metrics <- metrics %>% filter(method != 'SPADE') %>%
          mutate(method = fct_recode(method, 'DESpace2' = 'DESpace'))
metrics$method <- droplevels(metrics$method)

metrics <- metrics %>% dplyr::rename(RL = `relative length`,
                                     SS = `signal strength`) %>% 
  mutate(RL = factor(RL), 
         SS = factor(SS)) 

  
# metrics %>%
#   mutate(`RL` = scales::percent(`RL`))

p <- ggplot(metrics,
            aes(x = method, y = auprc, color = method)) +
  geom_jitter(width = 0.125, alpha = 0.7, size = 1.25) + 
  stat_summary(fun = mean, geom = "point", color = "black", size = 2, pch = 18) +
  # facet_grid(
  #   `signal strength` ~ `relative length`, 
  #   scales = "free_y",
  #   labeller = label_both
  # ) + 
  facet_wrap(
    ~ `SS` * `RL`,
    scales = "free_y",
    ncol = 3,
    labeller = \(x) label_both(x, multi_line = F, sep = ": ")
  ) +
  theme_minimal(base_family = "Arial", base_size = 6) +
  scale_color_manual(values = method_cols) +
  scale_x_discrete(
    # 1. Shrink outer padding around the first and last x categories
    expand = expansion(add = c(0.3, 0.3)),
    labels = function(x) str_wrap(x, width = 10)
  ) +
  labs(x = NULL, y = 'AUPRC') +
  theme(
    panel.grid = element_blank(),
    # Add border around each subpanel to create a grid between them
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    # 2. Controls panel width-to-height ratio: taller panels pull x-coordinates closer together
    aspect.ratio = 1,
    # 3. Tighten distance from plot border to x tick labels
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6, 
                               margin = margin(t = 2, unit = "pt")),
    axis.ticks.length.x = unit(2, "pt"),
    # 4. Tighten y-axis text/tick margins to minimize the gap between facet columns
    axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 6, 
                               margin = margin(r = .5, unit = "pt")),
    axis.ticks.length.y = unit(2, "pt"),
    axis.title.y = element_text(size = 6, margin = margin(r = 2, l = 0, unit = "pt")),
    legend.position = "none",
    strip.text = element_text(face = "plain", size = 6, margin = margin(t = 1, b = 1)),
    # 5. Inter-panel spacing
    panel.spacing.x = unit(5, "pt"),
    panel.spacing.y = unit(5, "pt"),
    plot.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "pt")
  )

print(p)
# Reduce output width from 5 to 4 to further compress x spacing across the 3 columns
ggsave(
  plot = p, 
  filename = glue('testing/figs/auprc_evaluation_w_spade.pdf'), 
  device = 'pdf', 
  width = 4.5, 
  height = 4.5, 
  units = 'in'
)



write_rds(metrics, 'testing/output/auprc_metrics_wo_spade.rds')
