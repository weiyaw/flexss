K <- list(pop = 4, sub = 4)
deg <- 1
size <- 1000
burn <- 0


## example with one population and multiple subject curves
## init <- list(pop = get_pls(data$y, Bmat$pop, Kmat) + rnorm(NCOL(Bmat$pop), sd = 100))
set.seed(1)
fm1 <- fit_bs_splines(simdata, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)
set.seed(1)
fm2 <- fit_bs_splines_v3(simdata, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)


## example with multiple population and subject curves
load_all()
set.seed(1)
fm3 <- fit_bs_splines(simdata2, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)

load_all()
set.seed(1)
fm4 <- fit_bs_splines_v3(simdata2, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)

