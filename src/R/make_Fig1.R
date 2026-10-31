source('src/argo/00-load-packages.R')

KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'

data_env <- read_rds('output/expression_data_env.rds')

data_env$cov_donor <- data_env$donor_identity
data_env$PLOT_POSTERIOR <- TRUE
data_env$PLOT_CONF_BANDS <- FALSE
data_env$RETURN_FITS <- FALSE
data_env$features_to_fit <- features_to_plot
data_env$kernel_params <- res.ARGO$params # if my parameters are already fit, I don't want to fit likelihood again

fit_ARGO(data_env)


features_to_plot <- funGP.df %>% arrange(desc(stats)) %>% 
  {bind_rows(slice_head(., n = 2), 
             slice_tail(., n = 2))} %>% pull(gene)

res.ARGO <- read_rds(glue('output/{CONTRAST}_N30_diff_testing_results.rds'))
funGP.df <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation), 
         fdr = p.adjust(exp(-neglogp), method = 'BH')
  ) %>%
  arrange(desc(neglogp)) %>%
  mutate(Rank = row_number())


require(readxl)
require(UpSetR)
library(PRROC)

res.list <- read_rds(glue('output/{CONTRAST}_DGE_results_list.rds'))
res.sig <- lapply(res.list, function(res) filter(res, fdr < 0.05) %>% pull(gene))
sapply(res.sig, length)

pdf('figs/upsetPlot_FDR_05.pdf', width = 4, height = 4)
theme_set(theme_gray(base_family = "Arial", 
                     base_size = 8))
print(upset(fromList(res.sig), 
            order.by = "freq", 
            nsets = length(res.sig),
            nintersects = 10, shade.alpha = .5, point.size = 2, line.size = .75,
            text.scale = c(
              .75, # 1. Intersection size bar numbers (top)
              .5, # 2. Intersection size axis title (Y-axis label)
              .25, # 3. Set size bar labels (left)
              .5, # 4. Set size axis title (X-axis label)
              .75, # 5. Set names (rows on the grid matrix)
              .75  # 6. Numbers above the dots (if turned on)
            )), 
      newpage= F)
dev.off()



res <- filter(funGP_res, abs(net_poe) > .25 & abs(net_poe) < .6) %>% 
  filter(fdr < .05) %>% arrange(desc(neglogp)) %>% 
  slice_head(., n=30)
# {bind_rows(slice_head(., n = 20), slice_tail(., n = 10))}

stat_lfc <- compute_stat_continuous(funGP_fits, N)
col_fun = colorRamp2(
  breaks = c(-1, 0, 1),                        # Min, Center, Max limits
  colors = c("#377eb8", "lightgrey", "#e41a1c")    # Blue to Red gradient
)

pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_partial_heatmap.pdf'), width = 4, height =3.5)
Heatmap(stat_lfc[res$gene, ], cluster_columns = F, name = 'LFC',
        show_row_names = T,
        show_column_names = F, 
        col = col_fun,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8), 
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 8, fontface = "bold"), # Title size 14 pt
          labels_gp = gpar(fontsize = 8)                    # Labels size 10 pt
        ))
dev.off()