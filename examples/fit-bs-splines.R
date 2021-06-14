if (!requireNamespace("magrittr", quietly = TRUE)) {
  stop("Package \"magrittr\" needed to run the examples. Please install it.",
       call. = FALSE)
}

library(magrittr)

## An example with multiple population, subject curves, and fixed effects
## load_all()
rm(list = ls())

## An example on simdata3. simdata3 is a simulated dataset from quadratic
## splines with 3 equidistant knots.
fit_data <- get_simdata3() %>%
  dplyr::mutate(pop = as.factor(pop),
                sub = as.factor(sub))

spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
            sub = ~s(x, by = sub, knots = 3, deg = 2,
                     intercept = TRUE, is_sub = TRUE, block_dim = 3))
fix <- y ~ pop + effect1 + effect2

set.seed(1)
fm3 <- fit_bsm(fixed = fix, data = fit_data, spline = spl,
                         size = 100, burn = 50, ridge = FALSE, init = NULL)

## predict works
predict(fm3)

## if you only need the population curves, pass in a 'level' argument to tell
## the function to only use the spline 'pop' and fixed/random effects.
predict(fm3, level = 'pop')

## if you want to predict at a new data point
newdata <- tibble::tribble(
  ~x, ~pop, ~sub, ~effect1, ~effect2,
  1, 1, 1, 1, 2,
  2, 2, 2, 0, 1,
  3, 1, 3, 0, 0,
  )
predict(fm3, dplyr::mutate_at(newdata, -1, as.factor))

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


## An example on the sitka data
spl <- list(ozone = ~s(days, knots = 4, deg = 2),
            id.num = ~s(days, by = id.num, knots = 4, deg = 2,
                        intercept = T, is_sub = T, block_dim = 3))

fm6 <- fit_bsm(fixed = log.size ~ 1, data = sitka, spline = spl,
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


## An example on the Berkeley growth data
growth <- fda::growth[c('hgtm', 'hgtf')] %>%
  purrr::map(~tibble::as_tibble(.x, rownames = 'age')) %>%
  purrr::map_dfr(~tidyr::pivot_longer(.x, -age, names_to = 'child', values_to = 'height'), .id = 'sex') %>%
  dplyr::mutate(age = as.numeric(age))

spl <- list(sex = ~s(age, by = sex, knots = 7, deg = 2),
            child = ~s(age, by = child, knots = 8, deg = 2, is_sub = T))
fm7 <- fit_bsm(height ~ sex, growth, spl, size = 100)

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







