## make AUPRC curves with genesets of interest
library(readxl)
library(PRROC)
KERNEL <- 'RBF'
CONTRAST <- 'old_vs_adult'
## load ranked lists from different DGE testing approaches
## compare with fungp results
df1 <- vroom(glue('data/processed/{KERNEL}_{CONTRAST}_funGP_testing.csv')) 
df2 <- vroom(glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv')) %>% arrange(desc(stats))
df3 <- vroom(glue('data/processed/DESeq_DGE_bin_{CONTRAST}.csv')) %>% arrange(desc(stats))

## read genesets
#-------------------------------------------------------------------------
## aging-signature genes from https://www.cell.com/cell/fulltext/S0092-8674%2823%2900805-X
#-------------------------------------------------------------------------
CAS <- vroom('data/sheets/cas-genes.txt', delim = '\t')
#----------------------------------------------------------------------------------------
## list of age-related differential genes from https://www.nature.com/articles/s41586-024-08350-8#Sec8
#----------------------------------------------------------------------------------------
## filtered list of genes with only restricted celltypes
neuronal.dge <- read_xlsx('data/sheets/cortex_striatum_age_DE_filtered.xlsx', sheet = 1)
glial.dge <- read_xlsx('data/sheets/cortex_striatum_age_DE_filtered.xlsx', sheet = 2)
names(neuronal.dge)[7] <- names(glial.dge)[7] <- 'label'
age.dge <- rbind(neuronal.dge, glial.dge)

# restrict to MERFISH panels
age.dge <- filter(age.dge, gene %in% df1$gene)
## filter out genes with dual trends across cell types
gene_summary <- age.dge %>% group_by(gene) %>%
  summarize(upreg = sum(age_effect_size > 0),
            downreg = sum(age_effect_size < 0),
            mixed = any(age_effect_size > 0) & any(age_effect_size < 0)
  ) %>% arrange(desc(mixed))

sum(gene_summary$mixed)
sum(gene_summary$upreg > 1)
sum(gene_summary$downreg > 1)

filter(gene_summary, gene %in% CAS$`Gene Symbol`)

## can filter out the mixed-signal genes if necessary
age.up <- filter(age.dge, (age_effect_size > 0) & 
                   (gene %in% (filter(gene_summary, upreg > 0)  %>% pull(gene))) )  %>% 
  distinct(gene, .keep_all = T)

age.dn <- filter(age.dge, (age_effect_size < 0) &
                   (gene %in% (filter(gene_summary, downreg > 0)  %>% pull(gene))) ) %>% 
  distinct(gene, .keep_all = T)

#------------------------------------------------------
## make the auprc curves
#------------------------------------------------------
compute_ranked_auprc <- function(ranked_list, signature){
  binary_labels <- ranked_list %in% signature
  scores <- rev(seq(length(ranked_list), 1))
  pg_score <- scores[binary_labels]
  ng_score <- scores[!binary_labels]
  
  pr_obj <- pr.curve(scores.class0 = pg_score, scores.class1 = ng_score, curve = T)
  return(pr_obj)
}


signature <- age.dn$gene
# signature <- CAS$`Gene Symbol`
## remove neg-shfit genes
# signature <- filter(deseq_dge, (gene %in% signature) & (stats < 0)) %>% pull(gene)

SX <- 'Jin et al aging down'
pr_fungp <- compute_ranked_auprc(df1 %>% arrange(desc(stats)) %>% pull(gene), signature)
pr_deseq_tissue <- compute_ranked_auprc(df2 %>% arrange(desc(stats)) %>% pull(gene), signature)
pr_deseq_bin <- compute_ranked_auprc(df3 %>% arrange(desc(stats)) %>% pull(gene), signature)

pdf(glue('figs/{CONTRAST}_AUPRC_domain_{SX}.pdf'), width = 6, height = 6)
par(mfrow = c(1,1), mar = c(8, 6, 4, 2) + 0.1, family = 'Arial', ps=16)
plot(pr_fungp, col = adjustcolor("red3", alpha.f = 1), main = glue("{SX} PR Curve"), 
     auc.main = F, lty=2)
plot(pr_deseq_tissue, color = adjustcolor("dodgerblue3", alpha.f = .5), add = TRUE)
plot(pr_deseq_bin, color = adjustcolor("forestgreen", alpha.f = 0.5), add = TRUE)

legend("topright", cex=.8,
       legend = c(paste("FunGP (AUPRC =", round(pr_fungp$auc.integral, 3), ")"),
                  paste("DESeq2 tissue (AUPRC =", round(pr_deseq_tissue$auc.integral, 3), ")"),
                  paste("DESeq2 bin (AUPRC =", round(pr_deseq_bin$auc.integral, 3), ")")),
       col = c("red3","dodgerblue3", "forestgreen"), 
       lty = c(2,1,1),             # Solid line type
       lwd = 2,             # Line width
       bty = "n")           # Removes the box border around the legend (optional)
dev.off()

#------------------------------------------------------
## identify genes uniquely identified by funGP
#------------------------------------------------------
df1 <- df1 %>% mutate(deseq_tissue_neglogp = df2$neglogp_combined,
                      deseq_tissue_minneglog = df2$neglogp_min,
                      deseq_bin_neglogp = df3$neglogp_combined)

df1sx <- filter(df1, ((neglogp > -log(.05)) & (abs(net_auec) > 0)) & 
                  
                  !(
                  # (deseq_tissue_neglogp > -log(.05)) |
                  (deseq_tissue_minneglog > -log(.05))
                  # (deseq_bin_neglogp > -log(.05))
                  )
                )

df1sx

max_val <- 2
myc <- colorRamp2(
  breaks = c(-max_val, 0, max_val), 
  colors = c("blue3", "lightgray", "red3")
)

pdf(glue('figs/{KERNEL}_{CONTRAST}_FC_heatmap_funGP_specific_vs_minlogp.pdf'), width = 6, height = 6)
Heatmap(stat_lfc[df1sx$gene, ], 
        show_column_names = F,
        cluster_columns = F, name = 'FC \nstatistic',
        col = myc,
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=4))
dev.off()

