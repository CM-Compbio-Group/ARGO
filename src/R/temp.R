# as_tibble(prop) %>% 
#   mutate(condition = bin_meta$condition) %>% 
#   group_by(condition) %>% 
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

# par(mfrow = c(1,2))
# plot(prop[, 'Vasc']~ bin_meta$condition)
# plot(prop[, 'Neuron']~ bin_meta$condition)
# 
# plot(prop[, 'Vasc']~ bin_meta$total_cells, pch = 20, col = factor(bin_meta$condition))
# plot(prop[, 'Neuron']~ bin_meta$total_cells, pch = 20, col = factor(bin_meta$condition))
# 
# plot(prop[, 'Neuron']~ bin_meta$LE, pch = 20, col = factor(bin_meta$condition))
# plot(prop[, 'Vasc']~ bin_meta$CTmvp, pch = 20, col = factor(bin_meta$condition))
# 
# 
# plot(prop[, 'cMESHyp']~ bin_meta$condition)
# plot(prop[, 'cMESInf']~ bin_meta$condition)
# plot(prop[, 'cMESVasc']~ bin_meta$condition)
# plot(prop[, 'cProlifMetab']~ bin_meta$condition)
# plot(prop[, 'cOXPHOS']~ bin_meta$condition)
# plot(prop[, 'Oligo']~ bin_meta$condition)
# plot(prop[, 'cNPC']~ bin_meta$condition)
# plot(prop[, 'cOPC1']~ bin_meta$condition)
# par(mfrow=c(1,1))

# # Variance-stabilizing transformation for multinomial proportions, centered per state
# prop <- prop[as.character(bin_meta$condition_key), ]
# prop_vst <- asin(sqrt(prop))
# prop_centered <- scale(prop_vst, center = TRUE, scale = FALSE)

# sample-1 is condition-1; sample-2 is condition-2
# bin_meta <- bin_meta %>% mutate(sample_id = fct_recode(condition, 
#                                                        'sample_1' = 'GBM',
#                                                        'sample_2' = 'IDHmutant'), 
#                                 condition_key =  condition_key)

# xx <- brain@meta.data %>% dplyr::select(all_of(ivy))
## clamp the enrichment scores at zero
# xx <- pmax(as.matrix(brain@meta.data %>% dplyr::select(all_of(ivy))), 0)


# par(mfrow = c(2, 2))
# for (col in colnames(scores_softplus)) {
#   boxplot(scores_softplus[, col] ~ brain$ivy_niche, 
#           main = paste(col, "(Softplus)"),
#           xlab = "Assigned Niche",
#           ylab = col,
#           outline = FALSE)
# }
# par(mfrow = c(1,1))

## bin the data 
# brain$bins <- bin_kdtree_adaptive(scores_softplus, target_spot_count = 500)
# pheatmap(prop, cluster_rows = T)

# I don't need to fit the model within each slide separately.
# I can just jointly model across all slides and split by condition -- not modeling donor effects


bin_cols <- colorRampPalette(brewer.pal(9, "Set3"))(N)
plotly::plot_ly(data = df, 
                x = ~CT, 
                y = ~CTpan, 
                z = ~LE, 
                color = ~bin, 
                colors = bin_cols,
                type = "scatter3d", 
                marker = list(size = 2.5, alpha=.5),
                mode = "markers")



## filter out bins that are extremely skewed b/w conditions
bin_tab <- table(df$bin, df$condition)
prop_bins <- bin_tab[, "GBM"] / rowSums(bin_tab)
hist(prop_bins)
sum((prop_bins < .05) | (prop_bins > .95)) #21 bins dominated by idh-wt spots

# filter out bins with extreme skew b/w conditions
# keep_bins <- which((prop_bins > .05) & (prop_bins < .95))
# df <- filter(df, bin %in% keep_bins)
# df$bin <- droplevels(df$bin)
# (N <- length(levels(df$bin)))

## now filter out bins with very few cells
hist(bin_meta$total_cells)
quantile(bin_meta$total_cells, c(.1,.2, .25,.5))
sum(bin_meta$total_cells < 20)

# drop bins with fewer than 50 spots
bin_meta <- bin_meta %>% filter(total_cells > 30)

## remove dropped bins
bin_meta$bin <- droplevels(bin_meta$bin)
(N <- length(levels(bin_meta$bin)))