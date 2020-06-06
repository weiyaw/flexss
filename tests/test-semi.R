# this script containts files for testing routines for fitting semiparametric
## models with some dummy design matrices.

library(magrittr)
load_all()

## test init_delta
rm(list = ls())
theta <- 1:3
init_theta(theta, 'a', 3, 1, NULL) # throw an error
init_theta(as.matrix(theta), 'a', 3, 1, NULL) # return theta with colnames
init_theta(as.matrix(theta), 'a', 3, 2, NULL) # throw an error
init_theta(as.matrix(theta), 'a', 2, 1, NULL) # throw an error
theta <- cbind(1:3, 2:4)
init_theta(theta, c('a', 'b'), 3, 2, NULL) # return theta with colnames
colnames(theta) <- c('a', 'b')
init_theta(theta, c('b', 'a'), 3, 2, NULL) # return theta with columns reversed
pls <- 10:12
init_theta(NULL, NULL, 3, 2, NULL) # throw an error
init_theta(NULL, NULL, 3, 1, pls) # return one column of pls
init_theta(NULL, NULL, 3, 2, pls) # return two colums of pls
init_theta(NULL, 'a', 3, 1, pls) # return two colums of pls with colnames
init_theta(NULL, c('a', 'b'), 3, 2, pls) # return two colums of pls with colnames
init_theta(NULL, c('a', 'b'), 3, 1, pls) # throw an error


## test init_delta
rm(list = ls())
delta <- 1:3
init_delta(delta, 'a', 3, 1)            # throw an error
init_delta(as.matrix(delta), 'a', 3, 1) # return delta with colnames
init_delta(as.matrix(delta), 'a', 3, 2) # throw an error
init_delta(as.matrix(delta), 'a', 2, 1) # throw an error
delta <- cbind(1:3, 2:4)
init_delta(delta, c('a', 'b'), 3, 2) # return delta with colnames
colnames(delta) <- c('a', 'b')
init_delta(delta, c('b', 'a'), 3, 2) # return delta with columns reversed
init_delta(NULL, NULL, 3, 1)         # return one column of gaussian rv
init_delta(NULL, NULL, 3, 2)         # return two colums of gaussian rv
init_delta(NULL, 'a', 3, 1) # return two colums of gaussian rv with colnames
init_delta(NULL, c('a', 'b'), 3, 2) # return two colums of gaussian rv with colnames
init_delta(NULL, c('a', 'b'), 3, 1) # throw an error

## test init_beta
rm(list = ls())
beta <- 1:3
init_beta(beta, 'a', 3)                     # throw an error
init_beta(beta, letters[1:3], 3)            # return beta with rownames
init_beta(t(beta), letters[1:3], 3)         # throw an error
init_beta(as.matrix(beta), letters[1:3], 3) # return beta with rownames
init_beta(as.matrix(beta), letters[1:3], 2) # throw an error
init_beta(NULL, 'a', 3)                     # throw an error
init_beta(NULL, letters[1:3], 3)            # return gaussian rv with rownames


## test update_prec_theta()
rm(list = ls())
## return a matrix with the gamma value on the non-zero cov entries, and NaN otherwise
x <- matrix(1:20, 5, 4)
## penalise the last 3 rows
set.seed(1)
update_with_gamma(x[3:5, ], -0.5, 0)
set.seed(1)
Kmat <- cbind(0, 0, diag(3))
update_prec_theta(x, Kmat, list(a = -0.5, b = 0)) / Matrix::crossprod(Kmat)

## penalise the 1st-order difference between rows of x
set.seed(1)
update_with_gamma(matrix(1, 4, 4), -0.5, 0)
set.seed(1)
Kmat <- matrix(c(-1, 1, 0, 0, 0,
                 0, -1, 1, 0, 0,
                 0, 0, -1, 1, 0,
                 0, 0, 0, -1, 1), 4, 5, TRUE)
update_prec_theta(x, Kmat, list(a = -0.5, b = 0)) / crossprod(Kmat) 


## test update_prec_delta()
## might need to change prior (i.e. v of wishart) to get proper conditional posteriors
## return a block diagonal matrix
rm(list = ls())
x1 <- matrix(1:12, 2, 6)
x2 <- matrix(1:18, 3, 6)

## if dim_block = 0, return diag(gamma)
set.seed(1)
update_with_gamma(rbind(x1, x2), -0.5, 0)
set.seed(1)
update_prec_delta(rbind(x1, x2), 0, list(a = -0.5, b = 0))

## if 0 < dim_block = 2 < dim_delta, return diag(wishart, diag(gamma))
set.seed(1)
update_with_wishart(x1, -1, diag(2))
update_with_gamma(x2, -0.5, 0)
set.seed(1)
update_prec_delta(rbind(x1, x2), 2, list(a = -0.5, b = 0, v = -1, lambda = diag(2)))

## if dim_block = dim_delta = 5, return wishart
set.seed(1)
update_with_wishart(rbind(x1, x2), -1, diag(5))
set.seed(1)
update_prec_delta(rbind(x1, x2), 5, list(v = -1, lambda = diag(5)))




## test update_prec_eps()
## return the same values between the two routines
rm(list = ls())
set.seed(1)
update_with_gamma(1:3, -0.5, 0)
set.seed(1)
update_prec_eps(c(4, 6, 8), 3:5, list(a = -0.5, b = 0))


## test update_prec_beta()
## return a diag matrix with the 1,3,5th diag element a gamma random variate
rm(list = ls())
set.seed(1)
update_with_gamma(c(11, 13, 15), -0.5, 0)
set.seed(1)
update_prec_beta(11:15, c(1, 3, 5), list(a = -0.5, b = 0))
update_prec_beta(11:15, c(1, 3, 5), NULL)                 # throw an error
update_prec_beta(11:15, NULL, NULL)                       # return a zero matrix
update_prec_beta(11:15, NULL, list(a = -0.5, b = 0))      # return a zero matrix
update_prec_beta(NULL, NULL, NULL)                        # return NULL
update_prec_beta(NULL, c(1, 3, 5), list(a = -0.5, b = 0)) # return NULL


## fit the model with some dummy data (i.e. fit with the truth)
## generate some dummy data
rm(list = ls())
set.seed(1)
Bmat <- list(pop = matrix(sample(1:5, 500 * 4, TRUE), 500, 4),
             sub = matrix(sample(1:5, 500 * 4, TRUE), 500, 4))
Xmat <- matrix(sample(1:5, 500 * 2, TRUE), 500, 2)
## each subject name must be unique, no nesting! 2 population and 5 subjects, 2
## subs in pop 1, 3 subs in pop 2, 100 samples for each sub
grp <- list(pop = c(rep(1, 200), rep(2, 300)),
            sub = rep(1:5, each = 100))

theta <- setNames(1:2, 1:2) %>% purrr::map(~rnorm(4))   # N(0, diag(4))
delta <- setNames(1:5, 1:5) %>% purrr::map(~rnorm(4)) # N(0, diag(4))
beta <- rnorm(2) # N(0, diag(2))

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
error <- rnorm(50) # N(0, diag(50))
y <- pop_term + sub_term + beta_term + error

## new Bmat
Bmat$pop <- purrr::map(unique(grp$pop), ~`[<-`(Bmat$pop, grp$pop != .x, , 0)) %>%
  {do.call(cbind, .)}

## new Kmat
Kmat <- Matrix::bdiag(diag(4), diag(4))

## true precision matrix, 0 entry for fixed effects
prec <- list(theta = Matrix::bdiag(diag(c(0, 0, 1, 1)), diag(c(0, 0, 1, 1))),
             delta = diag(4),
             beta = diag(c(0, 1)),
             eps = 1)

prior <- list(theta = list(a = -0.5, b = 0),
              delta = list(a = 0.5, b = 0, v = -1, lambda = diag(4)),
              beta = list(a = 0, b = 0),
              eps = list(a = -0.5, b = 0))

load_all()
set.seed(1)
## emperical bayes, fixed and random effects (ranef overriden by prec)
fm1 <- bayes_ridge_semi(y, grp, Bmat, Xmat, Kmat = Kmat, prec = prec, size = 200)
## full bayesian, fixed and random effects
fm2 <- bayes_ridge_semi(y, grp, Bmat, Xmat, Kmat = Kmat, dim_block = 0, ranef = c(1),
                        prior = prior, size = 200)
fm2$means
## full Bayesian, no fixed/random effects
load_all()
fm3 <- bayes_ridge_semi(y, grp, Bmat, NULL, Kmat = Kmat, 0, c(1), prior = prior, size = 200)
fm3$means
## emperical Bayesian, no fixed/random effects
fm4 <- bayes_ridge_semi(y, grp, Bmat, NULL, Kmat = Kmat, prec = prec, size = 200)



