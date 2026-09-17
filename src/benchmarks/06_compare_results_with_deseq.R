KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'

source('src/scripts/00-load-packages.R')
library(babelgene)
library(enrichR)

#--------------------------------------------
## plot deseq results
#--------------------------------------------
## compare with fungp results
df1 <- vroom(glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv')) 
df2 <- vroom(glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv')) %>% arrange(desc(stats))


goi <- c(setdiff(df1 %>% arrange(desc(stats)) %>% pull(gene) %>% head(n = 50), 
                 df2 %>% arrange(desc(stats)) %>% pull(gene) %>% head(n = 50)), 
         setdiff(df1 %>% arrange((stats)) %>% pull(gene) %>% head(n = 50), 
                 df2 %>% arrange((stats)) %>% pull(gene) %>% head(n = 50))
)

## arrange df2 in same order as df1
df2 <- df2 %>% arrange(match(gene, df1$gene))
df1 <- bind_cols(df1, df2 %>% select(neglogp_combined, neglogp_min, maxlogFC))

df1 <- df1 %>% mutate(status = case_when(
  deviation > 0 & maxlogFC > 0 ~ 'Upreg',
  deviation < 0 & maxlogFC < 0 ~ 'Downreg',
  TRUE ~ 'Discordant'
)) 


df1$to_label <- FALSE
df1$to_label[which(df1$gene %in% goi)] <- TRUE

#--------------------------------------------
## prepare volcano and label discordant genes by rank
#--------------------------------------------
pdf(glue('figs/{CONTRAST}_differential_significance_volcano.pdf'), width = 8, height = 6)
ggplot(df1, aes(x = maxlogFC, y = deviation)) +
  # 1. Plot points with grey color and 25% alpha transparency
  geom_point(aes(color = status), alpha = 0.35, size = 2, shape = 20) +
  scale_color_manual(
    values = c(
      "Upreg" = "#e41a1c",   # Professional Red
      "Downreg" = "#377eb8",  # Professional Blue
      "Discordant" = "gray60"     # Muted Gray
    )
  ) +
  geom_text_repel(
    data = . %>% filter(to_label),
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
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = .75) +
  # 3. Vertical effect size threshold lines
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = .75) +
  # 4. Axis Labels
  labs(
    x = "DESEq logFC ",
    y = "funGP max deviation (LFC)"
  ) +
  # xlim(c(-3.5, 3.5)) + 
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

#--------------------------------------------
## run enrichment tests w/ four sets of genes
#--------------------------------------------



#--------------------------------------------
#--------------------------------------------
df1 %>% 
  with(table(neglogp_combined > -log(.05), (neglogp > -log(.05) & abs(net_auec) > 0)))

## unique genes found by fungp
df1 %>% filter((neglogp > -log(.05) & abs(net_auec) > 0.25) & neglogp_combined < -log(.05)) %>%
  with(table(net_auec > 0))


## genes found by deseq
dfx <- df1 %>% filter(neglogp < -log(.05) & neglogp_min > -log(.05)) %>%
  arrange(desc(maxlogFC))
   # with(table(maxlogFC > 0))
  # %>% write_csv(file = 'data/processed/deseq_unique_dge.csv')
write_csv(dfx, 'data/processed/deseq_domain_minlogp_unique_dge.csv')

table(dfx$maxlogFC > 0)
## plot the fitted curves for these genes
genes_to_plot <- dfx$gene
source('src/scripts/0B-load-data.R')
kernelParams <- read_rds(glue('output/kernel_params_{KERNEL}.rds'))
params <- do.call('rbind', kernelParams)
colnames(params) <- c('bandwidth', 'sigma_f', 'sigma_d', 'sigma_n1', 'sigma_n2', 'sigma_n3')
FIT_DONOR_SWITCH = FALSE

## want to argue that genes identified w/ min-logp don't show strong signals
source('src/scripts/02_visualize_fitted_GPs.R')
ggsave( wrap_plots(plots[1:15], ncol = 4) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_posterior_fit_minlog_dfx_upreg.pdf')),
        width=18,height=15)

ggsave( wrap_plots(plots[16:21], ncol = 3) + 
          plot_layout(guides = "collect"), 
        filename = file.path('figs', glue('{KERNEL}_{CONTRAST}_posterior_fit_minlog_dfx_downreg.pdf')),
        width=12,height=8)

## visualize the lfc values from deseq for these genes
library(extrafont)
font_import(pattern = "Arial")
loadfonts(device = "pdf", quiet = TRUE)
max_val <- 2
my_colors <- colorRamp2(
  breaks = c(-max_val, 0, max_val), 
  colors = c("blue", "white", "red")
)

pdf(glue('figs/{CONTRAST}_DEseq_domain_minlogp_dx_LFC_heatmap.pdf'), width = 6, height = 6)
Heatmap(logFC[dfx$gene, ], name = 'LFC', 
        show_row_names = T,
        cluster_columns = F, cluster_rows = F, 
        col = my_colors,
        column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
        heatmap_legend_param = list(
          title = "LogFC",
          title_gp = gpar(fontsize = 16, fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
          legend_height = unit(4, "cm"), # Adjust height of the color bar
          grid_width = unit(5, "mm")      # Adjust width of the color bar
        ))
dev.off()


## how many significant genes by each test
df1 %>% filter((neglogp > -log(.05)) & (abs(net_auec) > .25)) %>% 
      with(table(net_auec > 0))

df2 %>% filter((neglogp_combined > -log(.05)) & (abs(maxlogFC) > 0)) %>% 
  with(table(maxlogFC > 0))

# df1 <- df1 %>% filter(neglogp > -log(.05) & (abs(net_auec) > .25)) 
             
# df1 %>% 
#   with(table(neglogp_combined > -log(.05)))

df1 %>% filter(neglogp > -log(.05)) %>%
  with(table(maxlogFC > 0, net_auec > 0))

table(df1$net_auec > 0, df1$maxlogFC > 0)

# df1 %>% filter((net_auec < 0) & (maxlogFC > 0)) 
# filter(df1, (maxlogFC > 0) & (net_auec <= 0)) %>% pull(gene)
# table(df1$net_auec > 0, df1$maxlogFC > 0)

#--------------------------------------------
#--------------------------------------------
## compare sign for genes
pdf(glue('figs/signed-rank-volcano_domain_{CONTRAST}.pdf'), width = 6, height = 6)
par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
plot(df2$maxlogFC~ df1$net_auec, ylab = 'DESeq domain-level LFC', 
     xlab = 'funGP net AUC', pch=19, cex=.75, 
     col=adjustcolor("gray3", alpha.f = 0.25),
     main = 'comparison of signed deviations')
abline(v = c(0,0), lty=2, lwd=1.5, col='red')
abline(h = c(0,0), lty=2, lwd=1.5,col='red')
dev.off()

#--------------------------------------------
## plot ranks for signature genes
#--------------------------------------------


#--------------------------------------------
## threshold fungp results and compare deseq calls for them
#--------------------------------------------
df2 <- filter(df2, neglogp_combined > -log(.05)) %>% arrange(desc(stats)) 
table(df2$maxlogFC > 0)

df1 <- df1 %>% arrange(desc(stats)) %>% {bind_rows(slice_head(., n = 102), slice_tail(., n = 94))}

df <- data.frame('gene' = c(df2$gene, df1$gene), 
                 'cluster' = rep(c('deseq_upreg', 'deseq_downreg', 
                                   'fungp_upreg', 'fungp_downreg'), c(102, 94, 102, 94)))
table(df$cluster)

df <- data.frame('gene' = c(df2 %>% arrange(desc(stats)) %>%
                              {bind_rows(slice_head(., n = 50), slice_tail(., n = 50))} %>%
                              pull(gene),
                            df1 %>% arrange(desc(stats)) %>%
                              {bind_rows(slice_head(., n = 50), slice_tail(., n = 50))} %>%
                              pull(gene)), 
                 'cluster' = rep(c('deseq_upreg', 'deseq_downreg', 
                                   'fungp_upreg', 'fungp_downreg'), each = 50))


dir.create(glue('output/enrichr/diffx/{CONTRAST}_domain_w_background'))

dbs <- c("GO_Biological_Process_2026", "TISSUES_Curated_2025", "MSigDB_Hallmark_2020",
         "CellMarker_2024", "GTEx_Aging_Signatures_2021", "Tabula_Muris",
         "Reactome_Pathways_2024")



for(clust in unique(df$cluster)){
  message(clust)
  genes <- filter(df, cluster == clust) %>% pull(gene)
  enr <- enrichr(genes, dbs)
  
  for(db in dbs){
    if(nrow(enr[[db]])){
      message(db)  
      pdf(glue('output/enrichr/diffx/{CONTRAST}_domain/{clust}_{db}.pdf'), width=10, height = 6)
      plot(plotEnrich(enr[[db]], showTerms = 5, numChar = 50, 
                      y  = "Count", 
                      # orderBy = "Adjusted.P.value", 
                      title = glue('{clust} {db}'))) + 
        theme(axis.text.y = element_text(size = 16)) # Change 12 to your desired font size
      
      dev.off()
    }
  }
}

#--------------------------------------------
## focus on genes missed by deseq but identified with fungp
#--------------------------------------------
# df1 <- filter(df1, abs(net_auec) > .25 & neglogp > -log(.05))
# df2 <- filter(df2, gene %in% df1$gene)
# all(df2$gene == df1$gene)
# table(df1$net_auec > 0, df2$maxlogFC > 0)
# plot(df1$Rank, df2$Rank)
# dfx <- df1[df2$neglogp_combined < -log(.05), ]


## uniquely identified genes
dfx <- filter(df1, (neglogp > -log(.05)) & (abs(net_auec) > .25) & 
                (neglogp_combined < -log(.05)))
table(dfx$net_auec > 0)

write_csv(dfx, glue('data/processed/{CONTRAST}_fungp_domain_unique_dge.csv'))

dfx <- dfx %>% arrange(desc(net_auec))

## visualize the heatmap for these genes
## set colorbar
# max_val <- max(abs(stat_lfc[dfx$gene, ]))
max_val <- 2
my_colors <- colorRamp2(
  breaks = c(-max_val, 0, max_val), 
  colors = c("blue", "white", "red")
)


pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_heatmap_funGP_dgx_domain.pdf'), width = 8, height = 8)
Heatmap(stat_lfc[dfx$gene, ], cluster_columns = F, cluster_rows = T,
        name = 'FC \nstatistic',
        col = my_colors,
        row_names_gp = gpar(fontsize=8), 
        column_names_gp = gpar(fontsize=4))
dev.off()

#-------------------------------------------
## check the deseq lfc values for these same genes
#-------------------------------------------
res.list <- read_rds('data/temp/90wk_vs_4wk_bin_level_deseq_res.rds')
p.vals <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(pvalue)))
logFC <- do.call('cbind', lapply(res.list, function(obs) obs %>% pull(log2FoldChange)))
rownames(logFC) <- rownames(p.vals) <- rownames(res.list[[1]])

max_val <- max(abs(logFC[dfx$gene, ]))
my_colors <- colorRamp2(
  breaks = c(-max_val, 0, max_val), 
  colors = c("blue", "white", "red")
)

pdf(glue('figs/{CONTRAST}_funGP_dgx_DEseq_domain_LFC_heatmap.pdf'), width = 8, height = 8)
Heatmap(logFC[dfx$gene, ], name = 'LFC', 
        show_row_names = T,
        cluster_columns = F, cluster_rows = T, 
        col = my_colors,
        column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
        heatmap_legend_param = list(
          title = "LogFC",
          title_gp = gpar(fontsize = 16, fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
          legend_height = unit(4, "cm"), # Adjust height of the color bar
          grid_width = unit(5, "mm")      # Adjust width of the color bar
        ))
dev.off()


pdf(glue('figs/{CONTRAST}_funGP_dgx_DEseq_domain_NEGLOGP_heatmap.pdf'), width = 8, height = 8)
Heatmap(-log(p.vals)[dfx$gene, ], name = 'neglogp', 
        show_row_names = T,
        cluster_columns = F, cluster_rows = T, 
        column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
        heatmap_legend_param = list(
          title = "LogFC",
          title_gp = gpar(fontsize = 16, fontfamily = "Arial"),
          labels_gp = gpar(fontsize = 10, fontfamily = "Arial"),
          legend_height = unit(4, "cm"), # Adjust height of the color bar
          grid_width = unit(5, "mm")      # Adjust width of the color bar
        ))
dev.off()

nlpcx <- apply(p.vals[dfx$gene, ], 1, cct_p.vals)
all(names(nlpcx) == dfx$gene)
plot(nlpcx, dfx$net_auec, pch=20, cex=1.25)


