y <- 1:16
grp <- list(pop = c(rep(1, 6), rep(2, 6), rep(3, 4)),
            sub = rep(1:8, each = 2))
Bmat <- list(pop = matrix(1:80, length(y), 5),
             sub = matrix(1:48, length(y),  3))
Kmat <- matrix(1:10, 2, 5)
dim_sub1 <- 3
burn <- 2
size <- 10
prior <- NULL

## some precalculation
rank_K <- NROW(Kmat)
n_terms_pop <- NCOL(Bmat$pop)
n_terms_sub <- NCOL(Bmat$sub)
n_samples <- NROW(Bmat$sub)
n_pops <- length(unique(grp$pop))
n_subs <- length(unique(grp$sub))

para <- list(Kmat = Kmat, dim_sub1 = dim_sub1, rank_K = rank_K,
             n_terms_pop = n_terms_pop, n_terms_sub = n_terms_sub,
             n_samples = n_samples, n_pops = n_pops, n_subs = n_subs)
subs_in_pop <- tapply(as.character(grp$sub), grp$pop, unique, simplify = FALSE)
idx_pop <- tapply(seq_len(n_samples), grp$pop, function(x) x, simplify = FALSE)
idx_sub <- tapply(seq_len(n_samples), grp$sub, function(x) x, simplify = FALSE)


## start of the test
## check check_grp, check_Bmat, check_prior
grp <- check_grp(grp)
Bmat <- check_Bmat(Bmat, Kmat)
prior <- check_prior(prior, dim_sub1)

## check calc_xB
all(calc_xB(Bmat, idx_sub, para)$pop[, , 2] == crossprod(Bmat$pop[3:4, ]))
all(calc_xB(Bmat, idx_sub, para)$subpop[, , 3] == crossprod(Bmat$sub[5:6, ], Bmat$pop[5:6, ]))
all(calc_xB(Bmat, idx_sub, para)$sub[, , 8] == crossprod(Bmat$sub[15:16, ]))

## check calc_Bxy
all(calc_Bxy(y, Bmat, idx_sub, para)$pop[, "3"] == crossprod(Bmat$pop[idx_sub[["3"]], ], y[idx_sub[["3"]]]))
all(calc_Bxy(y, Bmat, idx_sub, para)$sub[, "5"] == crossprod(Bmat$sub[idx_sub[["5"]], ], y[idx_sub[["5"]]]))

## check initialise_samples
initialise_samples(para, size, dim_sub1, grp)

## check initialise_with_pls
xKmat <- diag(ncol(Bmat$pop))
pls <- tcrossprod(solve(crossprod(Bmat$pop) + xKmat), Bmat$pop) %*% as.vector(y)
init <- initialise_with_pls(NULL, para, grp, pls) # produce init$pop with identical columns
initialise_with_pls(init, para, grp, NULL) # should give the same as the previous one






