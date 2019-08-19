K <- list(pop = 4, sub = 4)
deg <- 1
size <- 10
burn <- 5


## init <- list(pop = get_pls(data$y, Bmat$pop, Kmat) + rnorm(NCOL(Bmat$pop), sd = 100))
load_all()
fm <- fit_bs_splines(simdata, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)
