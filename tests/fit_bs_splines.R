## load_all()
rm(list = ls())

## An example with multiple population and subject curves. 'simdata2' is a
## simulated dataset from quadratic splines with 3 equidistant knots.
fit_data <- dplyr::mutate(get_simdata2(),
                          pop = as.factor(pop),
                          sub = as.factor(sub))

spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
            sub = ~s(x, by = sub, knots = 3, deg = 2,
                     intercept = TRUE, is_sub = TRUE, block_dim = 3))
fix <- y ~ pop

set.seed(1)
fm2 <- fit_bs_splines_v4(fixed = fix, data = fit_data, spline = spl,
                         size = 100, burn = 50, ridge = FALSE, init = NULL)

if (!all((fit_data$truth - predict(fm2)) < mean(abs(fit_data$truth)) * 0.1)) {
  warnings('Prediction varies too much from the truth.')
}

rm(list = ls())


## An example with multiple population, subject curves, and fixed
## effects. 'simdata3' is an extension of 'simdata2' with two extra fixed
## effects.
fit_data <- dplyr::mutate(get_simdata3(),
                          pop = as.factor(pop),
                          sub = as.factor(sub))

spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
            sub = ~s(x, by = sub, knots = 3, deg = 2,
                     intercept = TRUE, is_sub = TRUE, block_dim = 3))
fix <- y ~ pop + effect1 + effect2

set.seed(1)
fm3 <- fit_bs_splines_v4(fixed = fix, data = fit_data, spline = spl,
                         size = 100, burn = 50, ridge = FALSE, init = NULL)

if (!all((fit_data$truth - predict(fm3)) < mean(abs(fit_data$truth)) * 0.1)) {
  warnings('Prediction varies too much from the truth.')
}



