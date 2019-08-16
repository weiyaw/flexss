rm(list = ls())

## a random example (just to see whether the code runs)
y <- 1:16
grp <- list(pop = c(rep(1, 6), rep(2, 6), rep(3, 4)),
            sub = rep(1:8, each = 2))
Bmat <- list(pop = matrix(1:80, length(y), 5),
             sub = matrix(1:48, length(y),  3))
Kmat <- diag(ncol(Bmat$pop))
dim_sub1 <- 2
burn <- 2
size <- 10
prior <- NULL

bayes_ridge_sub_v2(y, grp, Bmat, Kmat, dim_sub1, burn, size, init = NULL, prior = NULL, verbose = TRUE)





