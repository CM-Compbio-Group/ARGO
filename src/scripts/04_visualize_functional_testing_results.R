#' **04_visualize_functional_testing_results** : use the funGP testing outputs to rank genes, and generate visualizations
#' a first-level look at funGP testing results
#' @params funGP_fits : list of fitted variables for each gene form funGP
#' @params x_observed, x_test : observed / test data points we fit GPs to
#' outputs : 1. volcano plots 
#'           2. ranked statistic - plot
#'           3. confidence bands for genes of interest
#'           4. posterior fits for genes of interest
#'           5. heatmaps indicating domains with significant deviations for each gene
#'           6. summary plots for exceedence statistics
#'           also writes summary of funGP testing as dataframe to file



library(EnhancedVolcano)
## visualize outputs of fungp tests
source('src/R/plot_conf_bands.R')
source('src/R/compute_exceedence_metrics.R')
source('src/R/compute_exceedence_statistics.R')

funGP_fits <- read_rds(glue('output/{KERNEL}_funGP_testing_{CONTRAST}_res.rds'))
load('data/temp/domain.rdata')

#----------------------------------------------------
## make deviation-heatmaps to summarize funGP findings
#----------------------------------------------------
stat_binary <- compute_stat_binary(funGP_fits, M)
pdf(glue('figs/{KERNEL}_{CONTRAST}_binary_heatmap.pdf'), width = 8, height = 8)
binary_colors <- c("1" = "red3","0" = "lightgrey", "-1" = "navy")
Heatmap(stat_binary, cluster_columns = F, name = 'binary \nstatistic',
        show_row_names = T,
        row_names_gp = gpar(fontsize=2),
        column_names_gp = gpar(fontsize=4),
        col = binary_colors,
        heatmap_legend_param = list(color_bar = "discrete")
)
dev.off()

#----------------------------------------------------

stat_lfc <- compute_stat_continuous(funGP_fits, M)
pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_heatmap.pdf'), width = 8, height =8)
Heatmap(stat_lfc, cluster_columns = F, name = 'FC \nstatistic',
        row_names_gp = gpar(fontsize=2),
        column_names_gp = gpar(fontsize=4))
dev.off()

#-------------------------------------
## gather results from tests
#-------------------------------------
## compute exceedence metrics
exc.metrics <- do.call('rbind', lapply(funGP_fits, compute_exceedence_metrics, x_test))

## aggregate results from functional testing
funGP_res <- funGP_fits %>%
  map_dfr(`[`, c('neglogp', 'deviation', 'index'), .id = 'gene') %>%
  bind_cols(exc.metrics) %>%
  mutate(across(c(neglogp, deviation, index), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         net_auec = auec_up - auec_down,
         net_poe = poe_up - poe_down) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 

#----------------------------------------------------
## summarize the exceedence results in a plot
pdf(glue('figs/{KERNEL}_{CONTRAST}_exceedence_stats_summary.pdf'), width = 8, height = 8)
par(mfrow = c(2,2),  mar = c(4, 4, 4, 4) + 0.25, family = 'Arial', ps=12)
with(funGP_res, plot(poe_up,poe_down, pch = 20, cex=1.5, 
                   col = adjustcolor(col='red', alpha=.5),
                   xlab = 'Percent domain up', 
                   ylab = 'Percent domain down'))

with(funGP_res, plot(neglogp~net_auec, pch = 20, cex=1.5, 
                   col = adjustcolor(col='red', alpha=.5),
                   ylab = 'neglogp analytic', 
                   xlab = 'net AUEC'))
abline(h = c(-log(.05), -log(.05)), lty=2, lwd=2, col='gray')
abline(v = c(-.25,-.25),  lty=2, lwd=2, col='gray')
abline(v = c(.25, .25),  lty=2, lwd=2, col='gray')

hist(funGP_res$net_auec, main = 'AUEC', xlab = '',
     col = adjustcolor(col='red', alpha=.5))

with(funGP_res, plot(net_auec~Rank, pch = 20, cex=1.5, 
                   col = adjustcolor(col='red', alpha=.5),
                   xlab = 'Gene Rank', 
                   ylab = 'net AUEC'))

dev.off()

#----------------------------------------------------
## DGE testing statistics at different effect-size threholds
thresh <- c(.25, .5, .75, 1)
summary_stats <- map_dfr(thresh, function(t) {
  # message(t)
  funGP_res %>% filter(neglogp > -log(.01) & abs(deviation) > t) %>%
    summarise(auec_threshold = t, 
              upregulated = sum(net_auec > 0),
              downregulated = sum(net_auec < 0),
              total = n())
})

summary_stats

## save to file
write_csv(funGP_res, file = glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv'))
#--------------------------------
## rank cell types by their statistic
#--------------------------------
pdf(glue('figs/{KERNEL}_{CONTRAST}_ranked_plot.pdf'), width = 6, height = 6)
ggplot(funGP_res %>% arrange(desc(stats)), aes(x = Rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
             color = "red", size = 2.5, alpha=0.5) +
  
  # Add the non-overlapping labels
  geom_text_repel(
    data = . %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
    aes(label = gene),
    size = 5, 
    family = 'Arial',
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "red", 
    max.overlaps = 50
  ) +
  # Labels and clean styling
  labs(
    x = "",
    y = "Signed neglogp",
    title = glue("funGP Rank Statistics")
  ) +
  theme_classic(base_family = 'Arial') +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size=14)
  )
dev.off()

#--------------------------------
## enhanced volcano plot
#--------------------------------
v_thresh <- 0.5
h_thresh <- -log(0.05)
pdf(glue('figs/{KERNEL}_{CONTRAST}_exceedence_volcano.pdf'), width = 8, height = 6)
funGP_res %>% mutate(status = case_when(
  neglogp > h_thresh & deviation > v_thresh ~ 'Upreg',
  neglogp > h_thresh & deviation < -v_thresh ~ 'Downreg',
  TRUE ~ 'Not Significant'
)) %>% arrange(status == 'Not Significant') %>%
ggplot(aes(x = deviation, y = neglogp)) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status), alpha = 0.35, size = 2, shape = 20) +
  scale_color_manual(
    values = c(
      "Upreg" = "#e41a1c",   # Professional Red
      "Downreg" = "#377eb8",  # Professional Blue
      "Not Significant" = "gray60"     # Muted Gray
    )
  ) + 
  geom_text_repel(
    data = . %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 10), slice_tail(., n = 10))},
    aes(label = gene),
    size = 5, 
    family = 'Arial',
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray", 
    max.overlaps = 20
  ) +
  # 2. Horizontal significance threshold line: -log10(0.05) ~ 1.301
  geom_hline(yintercept = h_thresh, color = "red", linetype = "dashed", size = .75) +
  # 3. Vertical effect size threshold lines
  geom_vline(xintercept = c(-v_thresh, v_thresh), color = "red", linetype = "dashed", size = .75) +
  # 4. Axis Labels
  labs(
    x = "max deviation (LFC)",
    y = expression(-log[10]~"analytic")
  ) + xlim(c(-3.5, 3.5)) + 
  # 5. Clean, professional theme minimalization
  theme_bw(base_family = 'Arial') +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size=14)
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1))
  )
dev.off()



#-------------------------------------
## visualize conf bands for significant genes
#-------------------------------------
top_diff <- funGP_res %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 5), 
                                  slice_tail(., n = 5)) %>% pull(gene)}

conf_band_plots <- lapply(top_diff, function(cell) 
                                    plot_CI_bands(cell, x = x_test, 
                                               band = funGP_fits[[cell]][['band']],
                                               stat = funGP_fits[[cell]][['statistic']], 
                                               alpha = 95)
)


ggsave( wrap_plots(conf_band_plots, ncol = 5) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_conf_bands.pdf')),
        width=18,height=8)

#-------------------------------------
## visualize posterior fits for these genes
#-------------------------------------
source('src/scripts/02_visualize_fitted_GPs.R')
plots <- visualize_fitted_GP(genes_to_plot = top_diff, kernel_params = params, FILTER_TO_CONTRAST = T)
ggsave( wrap_plots(plots, ncol = 5) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_posterior_fits.pdf')),
        width=18,height=8)
