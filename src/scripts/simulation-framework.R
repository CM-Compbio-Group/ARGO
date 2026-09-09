#' simulation framework for testing how good our algorithm is in detecting differential gene expression trends
#' simulate 2-condition scenario. $H_0 : u_1 = u_2, H_1 : \exists x \in X s.t. u_1(x) != u_2(x) $
#' execute this by drawing two latent functions $u_1, u_2 \in GP(0, \Sigma)$, and some time domain where we add introduce the deviation


## set dimensions of domain
T <- 100
x <- seq(0, 1, length.out = T)

## define parameters of covariance kernel
sigma_f <- .1
bw <- .075

kern_exponent <- dist(x)^2
