# compute extent of perturbation for each bandwidth

# perturb_amplitude <- c(0.3, 0.6, 1) # magnitude of LFC in perturbation
# perturb_bandwidth <- c(0.03, 0.06, 0.1) # width of the perturbation
# expt_df <- crossing(amplitude = perturb_amplitude,
#                     bandwidth = perturb_bandwidth,
#                     ) %>% 
#   mutate(Experiment = row_number())
dataset <- 'PANCREAS'
THRESH <- 0.2
expt_df$center <- runif(nrow(expt_df), 0.2, 0.8)
sce <- read_rds(glue('testing/data/raw/{dataset}_sce.rds'))

t_ <- sce$pseudotime
ord <- order(t_)

lfc_mat <- matrix(
  0,
  nrow = length(t_),
  ncol = nrow(expt_df),
)

for (j in seq_len(nrow(expt_df))) {
  lfc_mat[, j] <-
    expt_df$amplitude[j] *
    exp(
      -(t_ - expt_df$center[j])^2 /
        (2 * expt_df$bandwidth[j]^2)
    )
}

image(lfc_mat[ord, ])
lower_time_bound <- apply(lfc_mat, 2, function(z) min(t_[abs(z) > THRESH])) # min time point where absolute lfc exceeds THRESHold
upper_time_bound <- apply(lfc_mat, 2, function(z) max(t_[abs(z) > THRESH]) ) # upper bound of interval

time_bound_width <- upper_time_bound - lower_time_bound
plot(time_bound_width~ expt_df$bandwidth)
expt_df$time_bound_width <- time_bound_width

rel_length <- expt_df %>% group_by(bandwidth) %>% summarise(`relative length` = mean(time_bound_width))
