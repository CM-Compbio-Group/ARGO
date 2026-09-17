library(readxl)
library(PRROC)


# age.dge <- read_xlsx('data/sheets/age-DE_genes.xlsx', sheet='cluster table')

## keep only genes in MERFISH panel

#------------------------------------------------------
#------------------------------------------------------
genes_to_highlight <- CAS$`Gene Symbol`

pdf(glue('figs/{KERNEL}_{CONTRAST}_CAS_ranked_plot.pdf'), width = 6, height = 6.5)
ggplot(df2, aes(x = Rank, y = stats)) +
  # Draws the main sorted line/points
  geom_line(color = "grey70") +
  geom_point(color = "black", size = 1) +
  
  # Highlight the specific genes with a unique color
  geom_point(data = . %>% filter(gene %in% genes_to_highlight),
             color = "red", size = 2.5, alpha=0.5) +
  
  # Add the non-overlapping labels
  geom_text_repel(
    data = . %>% filter(gene %in% genes_to_highlight),
    aes(label = gene),
    size = 4,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "red", 
    max.overlaps = 50
  ) +
  
  # Labels and clean styling
  labs(
    x = "Rank",
    y = "Statistic Value",
    title = glue("{CONTRAST}\n funGP CAS")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12)
  )
dev.off()

#------------------------------------------------------
#------------------------------------------------------
filter(deseq_dge, gene %in% CAS$`Gene Symbol`)


