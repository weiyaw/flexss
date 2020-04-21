## this script containts files for testing routines for fitting semiparametric
## models with some dummy design matrices.

library(magrittr)
rm(list = ls())

## generate some dummy data (i.e. truth)
set.seed(1)
Bmat <- list(pop = matrix(sample(1:5, 500 * 4, TRUE), 500, 4),
             sub = matrix(sample(1:5, 500 * 4, TRUE), 500, 4))
Xmat <- matrix(sample(1:5, 500 * 2, TRUE), 500, 2)
## each subject name must be unique, no nesting! 2 population and 10 subjects, 4
## subs in pop 1, 6 subs in pop 2, 50 samples for each sub
grp <- list(pop = c(rep(1, 200), rep(2, 300)),
            sub = rep(1:10, each = 50))

theta <- setNames(1:2, 1:2) %>% purrr::map(~rnorm(4))   # N(0, diag(4))
delta <- setNames(1:10, 1:10) %>% purrr::map(~rnorm(4)) # N(0, diag(4))
beta <- rnorm(2) # N(0, diag(2))

## y = Bmat$pop * theta + Bmat$sub * delta + Xmat * beta + Gaussian error
## grp$pop and grp$sub must be sorted, o/w the resulting y will not match due to
## split.data.frame
ord <- order(grp$pop, grp$sub)
Bmat <- purrr::map(Bmat, ~.x[ord,])
Xmat <- Xmat[ord, ]
grp <- purrr::map(grp, ~.x[ord])

pop_term <- split.data.frame(Bmat$pop, grp$pop) %>%
  purrr::map2(theta, ~as.numeric(.x %*% .y)) %>%
  unlist(use.names = FALSE)
sub_term <- split.data.frame(Bmat$sub, grp$sub) %>%
  purrr::map2(delta, ~as.numeric(.x %*% .y)) %>%
  unlist(use.names = FALSE)
beta_term <- Xmat %*% beta
error <- rnorm(50) # N(0, diag(50))
y <- pop_term + sub_term + beta_term + error


## true precision matrix
prec <- list(pop = diag(4),
             sub = diag(4),
             beta = diag(2),
             eps = 1)

load_all()
sam <- bayes_ridge_semi(y, grp, Bmat, Xmat, prec = prec)
