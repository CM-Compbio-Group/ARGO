#============================================================
## find what resolution to fit data at by looking at gene expression sparsity levels
#============================================================
find_sparsity <- function(N, counts, metadata){
  bins <- ntile(metadata$pseudotime, N)
  slide_bin <- factor(paste(metadata$sample_id, bins, sep = '_'))
  bin_cts <- rowsum(counts, group = slide_bin)
  return(rowSums(bin_cts == 0))
}


n_range <- seq(10,90,10)
sparsity <- lapply(n_range, find_sparsity, counts, metadata)
sapply(sparsity, function(x) sum(x > 0))

par(mfrow = c(3,3))
for(i in seq_along(n_range)){
  hist(sparsity[[i]],  xlab = 'number of sparse genes per bin',
       main = glue('Sparsity with {n_range[i]} isodepth bins'))
}

sparsity_quantiles <- lapply(sparsity, function(x) quantile(x, c(.5, .75, .9, .99)))

t <- as.matrix(do.call('rbind', sparsity_quantiles))
rownames(t) <- n_range
print(t)
