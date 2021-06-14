if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package \"dplyr\" needed to run the test files. Please install it.",
       call. = FALSE)
}

if (!requireNamespace("withr", quietly = TRUE)) {
  stop("Package \"withr\" needed to run the test files. Please install it.",
       call. = FALSE)
}

## An example with multiple population and subject curves. 'simdata2' is a
## simulated dataset from quadratic splines with 3 equidistant knots.
test_that("Prediction varies from the truth (multi pop, no fixed effects)", {
  fit_data <- dplyr::mutate(get_simdata2(),
                            pop = as.factor(pop),
                            sub = as.factor(sub))
  
  spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
              sub = ~s(x, by = sub, knots = 3, deg = 2, is_sub = TRUE))
  fix <- y ~ pop
  
  withr::with_seed(10, {
    fm <- fit_bsm(fixed = fix, data = fit_data, spline = spl,
                            size = 100, burn = 50, ridge = FALSE, init = NULL)
  })
  
  expect_equal(fit_data$truth, predict(fm), tolerance = 0.03)
  expect_true(all((fit_data$truth - predict(fm)) < mean(abs(fit_data$truth)) * 0.1))
})



## An example with multiple population, subject curves, and fixed
## effects. 'simdata3' is an extension of 'simdata2' with two extra fixed
## effects.
test_that("Prediction varies from the truth (multi pop, with fixed effects)", {
  fit_data <- dplyr::mutate(get_simdata3(),
                            pop = factor(pop, levels = c('3', '2', '1')),
                            sub = as.factor(sub))
  
  spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
              sub = ~s(x, by = sub, knots = 3, deg = 2, is_sub = TRUE))
  fix <- y ~ pop + fixed1 + fixed2
  
  withr::with_seed(10, {
    fm <- fit_bsm(fixed = fix, data = fit_data, spline = spl,
                            size = 100, burn = 50, ridge = FALSE, init = NULL)
  })

  expect_equal(fit_data$truth, predict(fm), tolerance = 0.03) # as % of prediction
  expect_equal(predict(fm), fit_data$truth, tolerance = 0.03) # as % of truth
})


## An example with multiple population, subject curves, fixed and random
## effects. 'simdata4' is an extension of 'simdata2' with two extra fixed
## effects.
test_that("Prediction varies from the truth (multi pop, with fixed effects)", {
  fit_data <- dplyr::mutate(get_simdata4(),
                            pop = factor(pop, levels = c('3', '2', '1')),
                            sub = as.factor(sub))
  
  spl <- list(pop = ~s(x, by = pop, knots = 3, deg = 2),
              sub = ~s(x, by = sub, knots = 3, deg = 2, is_sub = TRUE))
  fix <- y ~ pop + fixed1
  
  withr::with_seed(10, {
    fm <- fit_bsm(fixed = fix, data = fit_data,
                            spline = spl, random = ~random1,
                            size = 100, burn = 50, ridge = FALSE, init = NULL)
  })
  
  expect_equal(fit_data$truth, predict(fm), tolerance = 0.04) # as % of prediction
  expect_equal(predict(fm), fit_data$truth, tolerance = 0.04) # as % of truth
})





