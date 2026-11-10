## plot simulated counts for a gene

mean_mu <- apply(mu_A, 2, mean)
hist(mean_mu)
q <- quantile(mean_mu, .85)
genes[mean_mu > q]

sort(mean_mu, decreasing = T)[1:10]

gene <- 'Tmsb4x'
plot(sce$pseudotime, mu_A[, gene], col = 'red3', pch = 20, ylim = range(c(mu_A[,gene], mu_B[,gene])))
points(sce$pseudotime, mu_B[, gene], col = 'blue3', pch = 20)

plot(sce$pseudotime, counts_A[[1]][gene, ], col = 'red3', pch = 20)
points(sce$pseudotime, counts_B[[1]][gene, ], col = 'blue3', pch=20)


df <- data.frame(
  pseudotime = rep(sce$pseudotime, 2),
  count      = c(counts_A[[1]][gene, ], counts_B[[1]][gene, ]),
  mu = c(mu_A[, gene], mu_B[, gene]),
  group      = rep(c("Group A", "Group B"), 
                   each = length(sce$pseudotime))
)

p <- ggplot(df, aes(x = pseudotime, color = group, fill = group)) +
  # Map raw counts to the scatter points
  geom_point(aes(y = count), alpha = 0.1, size = 0.8, shape = 16) +
  # Map mean parameters to the line (replaces geom_smooth)
  geom_line(aes(y = mu), linewidth = 0.8) +
  ylim(c(0, 70)) +
  scale_color_manual(values = c("Group A" = "red3", "Group B" = "blue3")) +
  scale_fill_manual(values = c("Group A" = "red3", "Group B" = "blue3")) +
  labs(
    x = "Pseudotime",
    y = paste0(gene, " Expression")
  ) + 
  theme_bw(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(2, "mm"),
    legend.box.spacing = unit(2, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

print(p)

ggsave(filename = glue('testing/figs/{gene}_scDesign_trend.pdf'), plot = p, width = 2, height = 2, units='in')

