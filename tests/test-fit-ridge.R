K <- list(pop = 4, sub = 4)
deg <- 1
size <- 1000
burn <- 0


## example with one population and multiple subject curves
## init <- list(pop = get_pls(data$y, Bmat$pop, Kmat) + rnorm(NCOL(Bmat$pop), sd = 100))
set.seed(1)
fm1 <- fit_bs_splines(simdata, K, deg, size, burn, ridge = FALSE, init = NULL, prior = NULL)

## rm(list = ls())
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
set.seed(1)
Xmat <- model.matrix(~effect1 + effect2, simdata3)[, -1] # col space of Xmat can't span 1
ranef <- NULL # all betas are fixed
load_all()
set.seed(1)
spl <- list(pop = ~s(x, by = pop, K = 4, deg = 1),
            sub = ~s(x, by = sub, K = 4, deg = 1, intercept = TRUE))
fix <- ~pop + effect1 + effect2

load_all()
set.seed(1)
fm4 <- fit_bs_splines_v3(simdata3, K = K, deg = deg,
                         spline = spl, fixed = fix,
                         ## Xmat = Xmat, ranef = ranef,
                         size = 100, burn = burn, ridge = FALSE,
                         init = NULL, prior = NULL)

load_all()
predict(fm4, data.frame(x = 1:3, pop = factor(1), sub = factor(1), effect1 = factor(1), effect2 = factor(2)))

plot(fm4$samples$coef$theta[2, 1, ])

predict_grp(diag(4), factor(c(1, 1, 2, 2)), matrix(1:8, 4, 2, dimnames = list(NULL, 1:2)))

plot_spline_v3(fm4)

