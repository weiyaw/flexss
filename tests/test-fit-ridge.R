K <- list(pop = 4, sub = 4)
deg <- 1
size <- 1000
burn <- 0


## example with one population and multiple subject curves
## init <- list(pop = get_pls(data$y, Bmat$pop, Kmat) + rnorm(NCOL(Bmat$pop), sd = 100))
set.seed(1)
fm1 <- fit_bs_splines(simdata, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)

load_all()
set.seed(1)
fm2 <- fit_bs_splines_v3(simdata, K = K, deg = deg,
                         size = size, burn = burn, ridge = FALSE, init = NULL, prior = NULL)
plot_spline_v3(fm2)


## example with multiple population and subject curves
load_all()
set.seed(1)
fm3 <- fit_bs_splines(simdata2, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)

load_all()
set.seed(1)
fm4 <- fit_bs_splines_v3(simdata2, K = K, deg = deg, Xmat = NULL, ranef = NULL,
                         size = size, burn = burn, ridge = FALSE, init = NULL, prior = NULL)
plot_spline(fm3)


## example with multiple population, subject curves, and fixed effects
load_all()
simdata3 <- get_simdata3()
set.seed(1)
Xmat <- model.matrix(~effect - 1, simdata3)
ranef <- NULL # all betas are fixed
fm4 <- fit_bs_splines_v3(simdata3, K = K, deg = deg, Xmat = Xmat, ranef = ranef,
                         size = size, burn = burn, ridge = FALSE, init = NULL, prior = NULL)
plot(fm4$samples$coef$theta[1, 1, ])
plot_spline(fm3)

