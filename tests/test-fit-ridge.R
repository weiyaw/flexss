library(magrittr)
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
load_all()
rm(list = ls())
set.seed(1)
spl <- list(pop = ~s(x, by = pop, knots = 4, deg = 2),
            sub = ~s(x, by = sub, knots = 4, deg = 2,
                     intercept = TRUE, is_sub = TRUE, block_dim = 3))
fix <- y ~ pop + effect1 + effect2
## prior <- list(spl = list(pop = list(a = -0.5, b = 0),
##                          sub = list(a = -0.5, b = 0,
##                                     v = -1, lambda = diag(2))),
##               eff = list(fixed = NULL),
##               eps = list(a = -0.5, b = 0))

## ## fit with nlme to get variance-cov matrix, and do emperical bayes to see whether the coef is correct
## Bmat <- get_design_bs(simdata3$x, K = 4, deg = 1)$design %*% get_transform_bs(4 + 1 + 1, 1 + 1)
## Xmat <- model.matrix(~ 0 + pop + effect1 + effect2, simdata3)
## bsub <- Bmat[, 1:2]
## isub <- Bmat[, -(1:2)]

## library(nlme)
## me1 <- lme(fixed = y ~ 0 + pop + effect1 + effect2 + bsub[, 2],
##            data = simdata3,
##            random = list(pop = pdIdent(~isub - 1),
##                          sub = pdBlocked(list(pdSymm(~bsub - 1),
##                                               pdIdent(~isub - 1)))))

## prec <- list(spl = list(pop = (as.matrix(me1$modelStruct$reStruct$pop) * me1$sigma^2) %>%
##                           {diag(c(0, 1 / diag(.)))},
##                         sub = solve(as.matrix(me1$modelStruct$reStruct$sub) * me1$sigma^2)),
##              eff = list(fixed = diag(0, 6)),
##              eps = 1 / me1$sigma^2)


## ## fit the (Bayesian) model with stan
## library(rstan)
## options(mc.cores = parallel::detectCores())
## rstan_options(auto_write = TRUE)
## stan_data <- list(J = length(unique(simdata3$sub)),
##                   N = length(simdata3$sub),
##                   ncol_X = NCOL(Xmat) + 1,
##                   ncol_Z = NCOL(Bmat) - 2,
##                   ncol_Xs = 2,
##                   ncol_Zs = NCOL(Bmat) - 2,
##                   y = simdata3$y,
##                   sub = as.numeric(as.character(simdata3$sub)),
##                   X = cbind(Xmat, Bmat[, 2]), Z = Bmat[, -(1:2)],
##                   Xs = Bmat[, 1:2], Zs = Bmat[, -(1:2)])

## bm1 <- stan(file = 'test-fit-ridge.stan',
##             data = stan_data,
##             iter = 2000, warmup = 1000, chains = 4,
##             control = list(adapt_delta = 0.8, max_treedepth = 16),
##             open_progress = F)

set.seed(1)
load_all()
fm5 <- fit_bs_splines_v4(fixed = fix, data = simdata3, spline = spl,
                         size = 4, burn = 0, ridge = FALSE, init = NULL)

load_all()
flatten_chain(fm5)

predict(fm5)
predict(fm5, level = 'pop')
newdata <- tibble::tribble(
  ~x, ~pop, ~sub, ~effect1, ~effect2,
  1, 1, 1, 1, 2,
  2, 2, 2, 0, 1,
  3, 1, 3, 0, 0,
)
predict(fm5, dplyr::mutate_at(newdata, -1, as.factor))

plot_dat <- simdata3 %>%
  tidyr::expand(x = seq(1, 20, length = 300),
                tidyr::nesting(pop, sub, effect1, effect2)) %>%
  dplyr::mutate(predpop = predict(fm5, newdata = ., level = 'pop'),
                predsub = predict(fm5, newdata = .))

ggplot(plot_dat) +
  geom_point(aes(x, y, col = pop), data = simdata3) +
  geom_line(aes(x, predsub, col = pop, group = sub)) +
  geom_line(aes(x, predpop, col = pop, group = sub), size = 1)

ggplot(plot_dat) +
  geom_point(aes(x, y, col = sub), data = simdata3) +
  geom_line(aes(x, predsub, col = sub, group = sub))


## sitka data
spl <- list(ozone = ~s(days, knots = 4, deg = 2),
            id.num = ~s(days, by = id.num, knots = 4, deg = 2,
                        intercept = T, is_sub = T, block_dim = 3))

fm6 <- fit_bs_splines_v4(fixed = log.size ~ 1, data = sitka, spline = spl,
                         size = 10)

plot_dat <- sitka %>%
  tidyr::expand(tidyr::nesting(id.num, ozone),
                days = seq(150, 674, length = 300)) %>%
  dplyr::mutate(predsub = predict(fm6, newdata = .),
                predpop = predict(fm6, newdata = ., level = 'ozone'),
                ozone = as.factor(ozone))

ggplot(dplyr::filter(plot_dat)) +
  geom_point(aes(days, log.size, col = factor(ozone)), data = sitka) +
  geom_line(aes(days, predsub, col = ozone, group = id.num)) +
  geom_line(aes(days, predpop, group = id.num), size = 2, col = 'black')


## growth data
growth <- fda::growth[c('hgtm', 'hgtf')] %>%
  purrr::map(~tibble::as_tibble(.x, rownames = 'age')) %>%
  purrr::map_dfr(~tidyr::pivot_longer(.x, -age, names_to = 'child', values_to = 'height'), .id = 'sex') %>%
  dplyr::mutate(age = as.numeric(age))

spl <- list(sex = ~s(age, by = sex, knots = 7, deg = 2),
            child = ~s(age, by = child, knots = 8, deg = 2, is_sub = T))
fm7 <- fit_bs_splines_v4(height ~ sex, growth, spl, size = 100)

plot_dat <- growth %>%
  tidyr::expand(tidyr::nesting(sex, child),
                age = seq(1, 18, length = 300)) %>%
  dplyr::mutate(predsub = predict(fm7, newdata = .),
                predpop = predict(fm7, newdata = ., level = 'sex'),
                sex = as.factor(sex))

ggplot(plot_dat) +
  geom_point(aes(age, height, col = child), data = growth) +
  geom_line(aes(age, predsub, col = child)) +
  geom_line(aes(age, predpop, group = sex), size = 3, col = 'black') +
  theme(legend.position = 'none') +
  facet_grid(cols = vars(sex))



