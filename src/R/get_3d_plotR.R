data_env <- read_rds('output/glioma_condition_data_env.rds')
# 2. Unpack everything from 'my_env' into the Global Environment
list2env(as.list(data_env), envir = .GlobalEnv)
all(metadata$grid_key == rownames(counts))
bin_meta <- cbind(metadata, counts)


res.ARGO <- read_rds('output/WT_vs_mutant_condition_results.rds')
res <- res.ARGO$resObjects %>%
  map_dfr(`[`, c('neglogp', 'deviation'), .id = 'gene') %>%
  mutate(across(c(neglogp, deviation), as.numeric)) %>% 
  mutate(stats = neglogp * sign(deviation)) %>%
  mutate(fdr = p.adjust(exp(-neglogp), method = 'BH')) %>%
  arrange(desc(neglogp)) %>% 
  mutate(Rank = row_number()) 

funGP_fits <- res.ARGO$resObjects

features <- res %>% arrange(stats) %>% {bind_rows(slice_head(., n = 5), slice_tail(., n = 5))} %>% pull(gene)

bin_WT <- filter(bin_meta, condition == 'GBM')
bin_Mutant <- filter(bin_meta, condition == 'IDHmutant')

plotly::plot_ly(data = bin_Mutant, 
                x = ~CTmvp, 
                y = ~CTpan, 
                z = ~LE, 
                color = as.formula(paste0("~", fx)), 
                # colors = bin_cols,
                type = "scatter3d", 
                marker = list(size =~sigDomainSize, 
                              opacity=0.15),
                mode = "markers")

for(fx in fetures){
  
  
  
  
  
}




# -------------------------------------------------------------------------
# 2. Build Plot 1: Mutant (Left Side)
# -------------------------------------------------------------------------

  # Significant points
   


# -------------------------------------------------------------------------
# 3. Build Plot 2: WT (Right Side)
# -------------------------------------------------------------------------
 

# -------------------------------------------------------------------------
# 4. Combine Side-by-Side and Configure Layout Domains
# -------------------------------------------------------------------------


## add a jitter to coordinates of IDH-mutant spots
ivy <- c('CT', 'CTmvp', 'CTpan', 'LE')
## find min distance b//w spots along each coordinate
coords <- bin_meta %>% select(all_of(ivy)) %>% distinct()
d <- quantile(dist(coords),.05)
jitter <- d 

rows_mutant <- which(bin_meta$condition == 'IDHmutant')
bin_meta[rows_mutant, ivy] <- bin_meta[rows_mutant, ivy] + jitter


states <- res$gene[1:5]

## every bin has unique coordinates -- has both mutant and WT specific data
# how do i visualize differences in 3D? 

## visualize abundance for cell types in 3d
bin_cols <- colorRampPalette(brewer.pal(9, "Set3"))(N)
