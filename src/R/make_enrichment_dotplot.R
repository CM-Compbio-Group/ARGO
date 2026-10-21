make_enrichment_dotplot <- function(gs, signature, cd1, cd2, wd = 4, ht = 6){
  
  gs$pathway <- str_replace_all(gs$pathway, pattern = '_', replacement = ' ')
  gs$pathway = str_wrap(gs$pathway, width = 30)
  
  ## keep only pathways of interest
  gs$pathway
  keep <- c(1,2,5,6,7,8,9,13,14)
  gs <- gs[keep, ]
  
  ## order pathways by pvals
  topPathwaysUp <- gs[ES > 0 & padj < .1][head(order(padj), n=20), pathway]
  topPathwaysDown <- gs[ES < 0 & padj < .1][head(order(padj), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  df <- data.frame(gs)
  rownames(df) <- df$pathway
  df <- df[topPathways, ]
  
  df$condition <- glue('{cd1} vs {cd2}')
  df$pathway <- factor(df$pathway, levels = rev(topPathways))
  
  pdf(glue('{signature}_{cd1}_vs_{cd2}.pdf'), width=wd, height = ht)
  plot(ggplot(data = df, aes(x = condition, y = pathway, group="",
                             color = `NES`, size = -log10(padj))) + 
         geom_point(alpha = 1) + 
         # ggtitle(signature) + 
         # scale_color_gradient2(low = "white",high = "red2", limits = (c(0, max(df$NES)))) +
         scale_color_gradient2(low = "blue",high = "red2", limits = range(df$NES)) +
         theme_bw() +
         theme(
           panel.grid = element_blank(),
           text = element_text(size = 16, family = 'Arial'),
           # axis.text.x = element_text(size=36),
           # axis.text.y = element_text(size=24),
           axis.title = element_text(size=24)) + 
         
         ylab("") + 
         xlab("") )
  
  dev.off()
  
}

make_enrichment_dotplot_rev <- function(gs, signature, cd1, cd2, wd = 4, ht = 6){
  
  gs$pathway <- str_replace_all(gs$pathway, pattern = '_', replacement = ' ')
  gs$pathway = str_wrap(gs$pathway, width = 30)
  
  ## order pathways by pvals
  topPathwaysUp <- gs[ES > 0 & padj < .05][head(order(padj), n=10), pathway]
  topPathwaysDown <- gs[ES < 0 & padj < .05][head(order(padj), n=10), pathway]
  topPathways <- c(topPathwaysDown, rev(topPathwaysUp))
  
  df <- data.frame(gs)
  rownames(df) <- df$pathway
  df <- df[topPathways, ]
  
  df$condition <- glue('Cluster {cd2} vs {cd1}')
  df$pathway <- factor(df$pathway, levels = rev(topPathways))
  df$NES <- -df$NES
  
  pdf(glue('{signature}_{cd2}_vs_{cd1}.pdf'), width=wd, height = ht)
  plot(ggplot(data = df, aes(x = condition, y = pathway, group="",
                             color = `NES`, size = -log10(padj))) + 
         geom_point(alpha = 1) + 
         # ggtitle(signature) +
         # scale_colour_manual(fill = 'blue3')+
         scale_color_gradient(low = "blue3",high = "red2") +
         theme(
           axis.text = element_text(size=16),
           axis.title = element_text(size=16)) + theme_bw() +
         ylab("") + 
         xlab("") )
  
  dev.off()
  
}

