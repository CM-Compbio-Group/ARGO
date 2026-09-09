# source('src/R/add_reflective_padding.R')
# load('data/temp/input_binned_data.rdata')
# cts <- assay(vsd)
# cts <- scale(t(assay(vsd)), center = T, scale = F)
# par(mfrow = c(1,1))
# boxplot(cts[,'C4b']~ bin_meta$age)
# boxplot(cts[,'B2m']~ bin_meta$age)
# boxplot(cts[,'Fth1']~ bin_meta$age)
# boxplot(cts[,'Il33']~ bin_meta$age)

# all(colnames(vsd) == bin_meta$slide_bin)
# genes <- rownames(vsd) # variable names
# ## append counts data to the metadata object
# bin_meta <- bin_meta %>% mutate(bin = as.character(bin)) %>% 
#   bind_cols(scale(t(assay(vsd)), center = T, scale = F)) 
# ## add padding to data
# bin_meta <- add_reflective_padding(bin_meta, n_bins = 8, max_bin = 25)
# plot(bin_meta$coord_x)
# write_rds(bin_meta, 'output/bin_meta.rds')
# sigx <- vroom('data/processed/sig_de_genes.csv')
# (genes_to_plot <- filter(sigx, shape_cluster==4) %>% pull(gene))
## do i want to extend test data? NO
# (delta <- median(diff(x_observed)))
# x_test <- c(-delta * rev(seq(n_bins)), x_test, 1+(delta*seq(n_bins)))
# M <- length(x_test)

#--------------------------------
## prepare projection matrix for observed data
## i've removed bins with fewer than 25 cells in a given slide
pp <- length(levels(bin_meta$bin))
qq <- length(unique(bin_meta$batch))
## get table of occurrences
tt <- c(table(bin_meta$bin, bin_meta$batch)) ## bins present per batch
## observed data is \hat{y} = P.y
P_obs <- matrix(0, nrow = nrow(bin_meta), ncol = pp*qq)
## fill out this matrix
ctx <- 1
for(ii in seq_along(tt)){
  if(tt[ii] == 1){
    P_obs[ctx, ii] <- 1
    ctx <- ctx +  1
  }
}

print(dim(P_obs))

#--------------------------------

## data is ordered as all slidexbins from age-4k, age-24wk and age-96wk respectively
#--------------------------------
## fit GPs to the gene type abundances
#--------------------------------
## no. of slides per condition
(slides_per_age <- bin_meta %>% distinct(batch, age) %>% dplyr::count(age) %>% pull(n))
n_slides <- sum(slides_per_age)

## prepare covariance structure for the donor-specific intercept 
donor_count_observed <- bin_meta %>%  count(donor_id, age) %>% arrange(age) %>% pull(n)
cov_donor <- Matrix::bdiag(lapply(donor_count_observed, function(x) matrix(1, nrow=x, ncol=x)))

donor_count <- bin_meta %>% distinct(batch, donor_id, age) %>% count(donor_id, age) %>% 
  arrange(age) %>% pull(n)


## define the co-regionalization matrix
W <- matrix(c(rep(c(1,0,0), slides_per_age[1]), 
              rep(c(0,1,0), slides_per_age[2]),
              rep(c(0,0,1),slides_per_age[3])), 
            nrow= n_slides, byrow = T)

## add scaling factor for correlation across slides per age group
rho <- 1
B <- rho * tcrossprod(W) + (1-rho) * diag(n_slides)

noise_diag <- c(rep(1, slides_per_age[1]),
                rep(2, slides_per_age[2]),
                rep(3, slides_per_age[3]))
