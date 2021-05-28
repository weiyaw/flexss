library(magrittr)

## Fit the model with some dummy design matrices (i.e. fit with the truth)

## Generate some dummy data
rm(list = ls())
set.seed(1)
Bmat <- list(pop = matrix(sample(1:5, 500 * 4, TRUE), 500, 4),
             sub = matrix(sample(1:5, 500 * 4, TRUE), 500, 4))
Xmat <- matrix(sample(1:5, 500 * 4, TRUE), 500, 4)
## each subject name must be unique, no nesting! 2 population and 5 subjects, 2
## subs in pop 1, 3 subs in pop 2, 100 samples for each sub
grp <- list(pop = c(rep(1, 200), rep(2, 300)),
            sub = rep(1:5, each = 100))


theta <- purrr::map(setNames(1:2, 1:2), ~rnorm(4)) # N(0, diag(4))
delta <- purrr::map(setNames(1:5, 1:5), ~rnorm(4)) # N(0, diag(4))
beta <- rnorm(4) # N(0, diag(4))

## y = Bmat$pop * theta + Bmat$sub * delta + Xmat * beta + Gaussian error
## grp$pop and grp$sub must be sorted, o/w the resulting y will not match due to
## split.data.frame
ord <- order(grp$pop, grp$sub)
Bmat <- purrr::map(Bmat, ~.x[ord,])
Xmat <- Xmat[ord, ]
grp <- purrr::map(grp, ~factor(.x[ord], unique(.x)))

pop_term <- split.data.frame(Bmat$pop, grp$pop) %>%
  purrr::map2(theta, ~as.numeric(.x %*% .y)) %>%
  unlist(use.names = FALSE)
sub_term <- split.data.frame(Bmat$sub, grp$sub) %>%
  purrr::map2(delta, ~as.numeric(.x %*% .y)) %>%
  unlist(use.names = FALSE)
beta_term <- Xmat %*% beta
truth <- pop_term + sub_term + beta_term
y <- truth + 0.1 * rnorm(length(truth)) # i.i.d. N(0, 0.1^2) error


## new Bmat
spl <- list()
spl$pop <- purrr::map(unique(grp$pop), ~`[<-`(Bmat$pop, grp$pop != .x, , 0)) %>%
  {do.call(cbind, .)}
attr(spl$pop, 'spl_dim') <- NCOL(spl$pop) / length(unique(grp$pop))
attr(spl$pop, 'is_sub') <- FALSE
attr(spl$pop, 'level') <- unique(grp$pop)
attr(spl$pop, 'penalty') <- diag(4)

spl$sub <- split.data.frame(Bmat$sub, grp$sub)
attr(spl$sub, 'spl_dim') <- NCOL(Bmat$sub)
attr(spl$sub, 'is_sub') <- TRUE
attr(spl$sub, 'level') <- as.character(unique(grp$sub))
attr(spl$sub, 'index') <- split(1:NROW(Bmat$sub), grp$sub)
attr(spl$sub, 'block_dim') <- 2

## eff1: random, eff2: fixed
eff <- list(eff1 = `colnames<-`(Xmat[, 1:2], c('a', 'b')),
            eff2 = `colnames<-`(Xmat[, 3:4], c('c', 'd')))

## true precision matrix, 0 entry for fixed effects
prec <- list(spl = list(pop = diag(c(0, 0, 1, 1)),
                        sub = diag(4)),
             eff = list(eff1 = diag(c(1, 1)),
                        eff2 = diag(c(0, 0))),
             eps = 1)

prior <- list(spl = list(pop = list(a = -0.5, b = 0),
                         sub = list(a = -0.5, b = 0,
                                    v = -1, lambda = diag(attr(spl$sub, 'block_dim')))),
              eff = list(eff1 = list(a = -0.5, b = 0),
                         eff2 = NULL),
              eps = list(a = -0.5, b = 0))


## A rough estimate of the magnitude of the true coeffcients. For testing if the
## posterior means is very different from the true coeffcients.
theta_scale <- mean(abs(do.call(cbind, theta)))
delta_scale <- mean(abs(do.call(cbind, delta)))
beta_scale <- mean(abs(beta))

fm1 <- bayes_ridge_semi_v4(y, spl, eff, prior = prior, prec = NULL, burn = 50, size = 100)
all(abs(do.call(cbind, theta) - fm1$means$coef$spline$pop) < theta_scale * 0.1)
all(abs(do.call(cbind, delta) - fm1$means$coef$spline$sub) < delta_scale * 0.1)
all(abs(beta - unlist(fm1$means$coef$effect)) < beta_scale * 0.1)

fm2 <- bayes_ridge_semi_v4(y, spl, eff, prior = NULL, prec = prec, burn = 50, size = 100)
all(abs(do.call(cbind, theta) - fm2$means$coef$spline$pop) < theta_scale * 0.1)
all(abs(do.call(cbind, delta) - fm2$means$coef$spline$sub) < delta_scale * 0.1)
all(abs(beta - unlist(fm2$means$coef$effect)) < beta_scale * 0.1)




