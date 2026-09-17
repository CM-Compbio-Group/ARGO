## find genes uniquely identified by fungp in domains
#-------------------------------------
## map fungp results to domamins
#-------------------------------------
df1 <- pvalsDf %>% filter(neglogp > -log(.05) & max_deviation > .5 )
test_genes <- df1$gene

## find sig rois for these genes and map into tissue domains
roi <- lapply(exact_pvals[test_genes], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
lx <- sapply(roi, length)
hist(lx)

## extract metadata from all slides
dfm <- bind_rows(lapply(dxlist, function(dx) dx@meta.data))
## convert to percent of spots in domain
# what percent of spots in the tissue domain are spanned by the roi?
n_spots_domain <- table(dfm$tissue)

## no. of spots in roi per tissue domain
n_spots_roi_per_domain <- do.call('rbind', lapply(roi, function(rx) {
  filter(dfm, binned_isodepth %in% rx) %>%
    with(table(tissue))
})
)

pct_roi_per_domain <- apply(n_spots_roi_per_domain, 1, function(x) x/n_spots_domain)
pheatmap(t(pct_roi_per_domain))


## visualize genes w/out global trends
goi <- filter(df1, net_poe < .6) %>% pull(gene)
pheatmap(t(pct_roi_per_domain)[goi, ])

## which of these are not significant by deseq?
filter(df2, gene %in% goi)
pheatmap(p.vals[goi, ] < .05)

#-------------------------------------
## consider significant if span is over 50% of domain?
#-------------------------------------
mask_fges <- t(pct_roi_per_domain) > .5
pheatmap(mask_fges)

colSums(mask_fges)
tx <- rowSums(mask_fges)
plot(lx ~ tx)

## there are genes spanning all eight domains w/ very low no. of rois
df1$n_domain <- tx
df1$n_roi <- lx

## genes spanning all domains with few rois
filter(df1, tx == 8 & lx < 70)

#-------------------------------------
## read deseq data
#-------------------------------------
df2 <- vroom(glue('data/processed/DESeq_DGE_tissue_{CONTRAST}.csv')) %>% arrange(desc(stats))
df2 <- df2 %>% mutate(neglogp = pmax(neglogp_combined, neglogp_min)) %>% 
  filter(neglogp > -log(.05) & abs(maxlogFC) > 0.25)

load(glue('output/DESeq_domain_{CONTRAST}_results.rdata'))
## make mask for deseq results
all(colnames(p.vals) == colnames(t(pct_roi_per_domain)))
mask_deseq <- p.vals[df1$gene, ] < 0.05 

#-------------------------------------
## find discordant pairs
#-------------------------------------
mask_discordant <- mask_fges - mask_deseq
mask_discordant <- mask_discordant[apply(abs(mask_discordant), 1, sum) > 0, ]

## keep only ones w/ positive shift
mask_discordant <- mask_discordant[rowSums(mask_discordant) > 0, ]
pheatmap(mask_discordant)

#-------------------------------------
## filter genes to visualize based on inter-donor variance
#-------------------------------------
df1 <- filter(df1, gene %in% rownames(mask_discordant))
## within each isodepth bin, i compute sd across all donors in age group
## take median sd across bins per age group 
return_age_sd <- function(gene){
  bin_meta %>% group_by(bin, age) %>%
    summarize(sd_val = sd(!!sym(gene), na.rm = TRUE), .groups = "drop")  %>%
    group_by(age) %>% summarize(sd_age = median(sd_val)) %>% pull(sd_age)
}
  
genes <- df1$gene
age_donor_sd <- do.call('rbind', map(genes, return_age_sd))
rownames(age_donor_sd) <- genes
colnames(age_donor_sd) <- c('sd_4wk','sd_24wk', 'sd_90wk')

df1 <- bind_cols(df1, age_donor_sd)
df1$donor_sd <- apply(age_donor_sd, 1, max)

## threshold genes by inter-donor variation
quantile(df1$donor_sd, c(.5,.75,.8, .9))

df1 %>% arrange(desc(deviation))

plot(df1$donor_sd~ df1$Rank)
genes_to_plot <- filter(df1, gene %in% rownames(mask_discordant) & 
                          donor_sd < quantile(donor_sd, .9)) %>% 
  arrange(desc(deviation)) %>%
  {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))} %>% pull(gene)

genes_to_plot

mask_discordant[genes_to_plot, ]
mask_deseq[genes_to_plot, ]
age_donor_sd[genes_to_plot, ]

#-------------------------------------
## visualize genes in space
#-------------------------------------
## assign regions from discordant mask
dmsk <- mask_discordant[genes_to_plot, ]
apply(dmsk, 1, which.max)

xt <- colnames(tt)[apply(dmsk, 1, function(t) order(t, decreasing=T)[1])]
names(xt) <- genes_to_plot


dfm <- filter(dfm, age %in% c('24wk', '90wk'))
## extract counts for genes of interest across all slides
cts <- do.call('cbind', lapply(dxlist[3:6], function(dx) LayerData(dx)[genes_to_plot, ]))
isodepth_bins <- do.call('c', lapply(dxlist[3:6], function(dx) dx[['binned_isodepth', drop=T]]))
domainx <- dfm%>% pull(tissue)
roi <- lapply(exact_pvals[genes_to_plot], function(kxx) x_test[which(abs(kxx$statistic) > kxx$band)])
sapply(roi, length)


## interesting genes found!!