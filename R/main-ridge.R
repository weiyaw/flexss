block_diag <- function(..., size = NULL) {

    ## Construct a big matrix with its diagonal elements the matrices provided
    ## in "A". If "A" is a matrix, the function returns a matrix with "size" of
    ## "A" in the diagonal elements. If "A" is a list, it returns a matrix with
    ## diagonal elements as all the matrices contained inside the list.

    A <- list(...)
    if (length(A) == 1 && is.matrix(A[[1]])) {
        if (is.null(size)) {
            stop("Size of the resulting matrix not supplied.")
        }
        n_row <- NROW(A[[1]])
        n_col <- NCOL(A[[1]])
        row_idx <- seq(1, n_row)
        col_idx <- seq(1, n_col)
        res <- matrix(0, size * n_row, size * n_col)
        for (i in seq(0, size - 1)) {
            res[i * n_row + row_idx, i * n_col + col_idx] <- A[[1]]
        }
        return(res)

   } else if (length(A) > 1) {
       if (!all(vapply(A, is.matrix, TRUE))) {
            stop("The list contains non-matrix objects.")
       }
        dims <- vapply(A, dim, c(1, 1))
        total_dims <- rowSums(dims)
        res <- matrix(0, total_dims[1], total_dims[2])
        row_rolling <- col_rolling <- 0
        for (i in seq(1, NCOL(dims))) {
            row_idx <- row_rolling + seq(1, dims[1, i])
            col_idx <- col_rolling + seq(1, dims[2, i])
            res[row_idx, col_idx] <- A[[i]]
            row_rolling <- row_rolling + dims[1, i]
            col_rolling <- col_rolling + dims[2, i]
        }
        return(res)

   } else {
       warning("Non-matrix or list object supplied")
   }
}

## source("subor.R")
## This is a Gibbs sampler for Bayesian ridge; see thesis.
## Requirements: Bmat, y, Kmat
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge <- function(y, Bmat, Kmat, burn, size, init = NULL, verbose = TRUE) {

    ## hyperparemeters for priors
    ig_a_pop <- 0.001                        # flat prior on standard deviations
    ig_b_pop <- 0.001
    ig_a_eps <- 0.001                        # flat prior on standard deviations
    ig_b_eps <- 0.001

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    xBmat <- crossprod(Bmat)
    xKmat <- crossprod(Kmat)
    Bxy <- crossprod(Bmat, y)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    precision = list())

    ## initialise theta using ols
    if (is.null(init)) {
        kcoef_pop <- tcrossprod(solve(crossprod(Bmat)), Bmat) %*% as.vector(y)
    } else {
        if (length(init$pop) == n_terms) {
            kcoef_pop <- as.vector(init$pop)
        } else {
            stop("Invalid initial values.")
        }
    }

    for (k in seq.int(-burn + 1, size)) {

        ## update sigma^2_theta
        shape_pop <- 0.5 * rank_K + ig_a_pop
        rate_pop <- 0.5 * crossprod(Kmat %*% kcoef_pop) + ig_b_pop
        kprec_pop <- rgamma(1, shape = shape_pop, rate = rate_pop)

        ## update sigma^2_epsilon
        shape_eps <- 0.5 * n_samples + ig_a_eps
        rate_eps <- 0.5 * crossprod(y - Bmat %*% kcoef_pop) + ig_b_eps
        kprec_eps <- rgamma(1, shape = shape_eps, rate = rate_eps)

        ## update theta
        M <- chol2inv(chol(xBmat + kprec_pop / kprec_eps * xKmat))
        mu <- M %*% Bxy
        sig <- 1 / kprec_eps * M
        kcoef_pop <- t(mvtnorm::rmvnorm(1, mu, sig))

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        if (k > 0) {
            samples$precision$pop[k] <- kprec_pop
            samples$precision$eps[k] <- kprec_eps
            samples$population[, k] <- kcoef_pop
        }
    }
    means <- list(population = rowMeans(samples$population))
    list(means = means, samples = samples)
}

## initialise theta with a penalised LS estimate, and delta with rnorms with small sd
initialise_with_pls <- function(init, n_terms, grp, pls) {
    ## init: NULL or list(pop, sub); any of the element can be NULL
    n_subs <- length(unique(grp))

    if (is.null(init$pop)) {
        init$pop <- pls
    } else {
        if (length(init$pop) == n_terms) {
            init$pop <- as.vector(init$pop)
            cat("Population initial values supplied.\n")
        } else {
            stop("Invalid dimension of population initial values.")
        }
    }
    if (is.null(init$sub)) {
        init$sub <- matrix(rnorm(n_terms * n_subs) * 0.01, n_terms, n_subs,
                           dimnames = list(NULL, levels(grp)))
    } else {
        if (dim(init$sub) == c(n_terms, n_subs)) {
            init$sub <- as.matrix(init$sub)
            cat("Subjects initial values supplied.\n")
            if (is.null(colnames(init$sub))) {
                colnames(init$sub) <- levels(grp)
            } else {
                init$sub <- init$sub[, levels(grp)]
            }
        } else {
            stop("Invalid dimension of subject initial values.")
        }
    }
    init
}

## calculates unnormalised log posterior for the subject specific model
## vector of coef_pop
## matrix of coef_sub
## list of prec: pop, sub1, sub2, eps
## list of contrib: pop, sub
## list of para: xKmat, rank_K, dim_sub1, y
## list of prior hyperparameters: ig_a, ig_b, iw_v, iw_lambda
logprior_sub <- function(coef_pop, coef_sub, prec, para, prior) {

    check_names <- prod(c("pop", "sub1", "sub2", "eps") %in% names(prec),
                        c("xKmat", "rank_K", "dim_sub1") %in% names(para),
                        c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))

    if (check_names == 0) {
        stop("Missing elements in the list.")
    }

    xKmat <- para$xKmat
    rank_K <- para$rank_K
    dim_sub1 <- para$dim_sub1
    ig_a <- prior$ig_a
    ig_b <- prior$ig_b
    iw_v <- prior$iw_v
    iw_lambda <- prior$iw_lambda          # assuming symmetric

    prec_sub <- block_diag(prec$sub1, diag(prec$sub2, NROW(coef_sub) - dim_sub1))

    ## inverse gamma
    lp_prec_pop <- (ig_a$pop + 1) * log(prec$pop) - ig_b$pop * prec$pop
    lp_prec_sub2 <- (ig_a$sub2 + 1) * log(prec$sub2) - ig_b$sub2 * prec$sub2
    lp_prec_eps <- (ig_a$eps + 1) * log(prec$eps) - ig_b$eps * prec$eps

    ## inverse Wishart
    lp_prec_sub1 <- 0.5 * (iw_v + dim_sub1 + 1) * determinant(prec$sub1)$modulus -
       0.5 * sum(diag(iw_lambda %*% prec$sub1))

    ## Gaussian
    lp_pop <- 0.5 * rank_K * log(prec$pop) -
        0.5 * prec$pop * crossprod(coef_pop, xKmat %*% coef_pop)
    lp_sub <- NCOL(coef_sub) / 2 * determinant(prec_sub)$modulus -
        0.5 * sum(apply(coef_sub, 2, function(x) crossprod(x, prec_sub %*% x)))

    lp_prec_pop + lp_prec_sub1 + lp_prec_sub2 + lp_prec_eps + lp_pop + lp_sub
}

## calculates unnormalised log-likelihood for the subject specific model
## list of prec: eps
## resids vector
## n_samples
loglike_sub <- function(prec, resids, n_samples) {
    check_names <- prod(c("eps") %in% names(prec))

    if (check_names == 0) {
        stop("Missing elements in the list.")
    }

    n_samples / 2 * log(prec$eps) - 0.5 * prec$eps * crossprod(resids)
}

## update precision
update_prec <- function(coef, resids, para, prior) {

    ## coef: list(pop = n_terms * 1 vec, sub = n_terms * n_subs mat)
    ## resid: n_samples * 1 vec
    ## para: list(Kmat, rank_K, dim_sub1, n_samples, n_subs, n_terms)
    ## prior: list(ig_a = list(pop, sub2, eps), ig_b = list(pop, sub2, eps),
    ##        iw_v, iw_lambda = dim_sub1 * dim_sub1 mat)

    Kmat <- para$Kmat
    rank_K <- para$rank_K
    n_samples <- para$n_samples
    n_subs <- para$n_subs
    n_terms <- para$n_terms
    dim_sub1 <- para$dim_sub1

    ig_a <- prior$ig_a                  # a list
    ig_b <- prior$ig_b                  # a list
    iw_v <- prior$iw_v                  # a number
    iw_lambda <- prior$iw_lambda        # a matrix, assuming symmetric

    prec <- list(pop = NA, eps = NA, sub1 = NA, sub2 = NA)

    ## update sigma^2_epsilon
    shape_eps <- 0.5 * n_samples + ig_a$eps
    rate_eps <- 0.5 * crossprod(resids) + ig_b$eps
    prec$eps <- rgamma(1, shape = shape_eps, rate = rate_eps)

    ## update sigma^2_theta
    shape_pop <- 0.5 * rank_K + ig_a$pop
    rate_pop <- 0.5 * crossprod(Kmat %*% coef$pop) + ig_b$pop
    prec$pop <- rgamma(1, shape = shape_pop, rate = rate_pop)

    ## update Sigma_dev1
    df_sub1 <- iw_v + n_subs
    scale_sub1 <- iw_lambda + tcrossprod(coef$sub[1:dim_sub1, , drop = FALSE])
    inv_scale_sub1 <- chol2inv(chol(scale_sub1))
    prec$sub1 <- rWishart(1, df = df_sub1, Sigma = inv_scale_sub1)[, , 1]

    ## update sigma^2_dev2
    shape_sub2 <- 0.5 * n_subs * (n_terms - dim_sub1) + ig_a$sub2
    rate_sub2 <- 0.5 * crossprod(c(coef$sub[-(1:dim_sub1), ])) + ig_b$sub2
    prec$sub2 <- rgamma(1, shape = shape_sub2, rate = rate_sub2)

    prec
}


## This is a Gibbs sampler for longitudinal Bayesian ridge; see thesis.
## Requirements: Bmat, y, grp, Kmat, dim_sub1
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_sub <- function(y, grp, Bmat, Kmat, dim_sub1, burn, size,
                            init = NULL, prior = NULL, verbose = TRUE) {

    if (is.null(prior)) {
        prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      iw_v = dim_sub1 + 1,
                      iw_lambda = diag(dim_sub1)) # assuming symmetry
    }

    check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))
    if (check_names == 0) {
        stop("Missing prior hyperparameters.")
    }

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    idx <- tapply(seq_len(n_samples), grp, function(x) x)

    ## some crossproducts precalculation
    xBmat <- crossprod(Bmat)
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
    }
    xKmat <- crossprod(Kmat)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)),
                    lp = rep(NA, size),
                    ll = rep(NA, size))

    ## initialise theta with a penalised LS estimate, and delta with rnorms with
    ## small sd
    pls <- tcrossprod(solve(crossprod(Bmat) + xKmat), Bmat) %*% as.vector(y)
    init <- initialise_with_pls(init, n_terms, grp, pls)
    kcoef <- list(pop = init$pop, sub = init$sub)

    ## initialise prediction contribution by population coefs and subjects deviations
    kcontrib_pop <- Bmat %*% kcoef$pop
    kcontrib_sub <- rep(NA, n_samples)
    for (i in levels(grp)) {
        kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
    }

    for (k in seq.int(-burn + 1, size)) {
        ## update precisions
        para <- list(Kmat = Kmat, rank_K = rank_K, dim_sub1 = dim_sub1,
                     n_samples = n_samples, n_subs = n_subs, n_terms = n_terms)
        kresids <- y - kcontrib_pop - kcontrib_sub
        kprec <- update_prec(kcoef, kresids, para, prior)

        ## update theta
        M_pop <- chol2inv(chol(xBmat + kprec$pop / kprec$eps * xKmat))
        mu <- M_pop %*% crossprod(Bmat, y - kcontrib_sub)
        sig <- 1 / kprec$eps * M_pop
        kcoef$pop <- t(mvtnorm::rmvnorm(1, mu, sig))
        kcontrib_pop <- Bmat %*% kcoef$pop

        ## update delta
        kprec_sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms - dim_sub1))
        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec$eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec$eps
            kcoef$sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec$pop
            samples$precision$sub1[, , k] <- kprec$sub1
            samples$precision$sub2[k] <- kprec$sub2
            samples$precision$eps[k] <- kprec$eps
            samples$population[, k] <- kcoef$pop
            samples$subjects[, , k] <- kcoef$sub

            ## calculate unnormalised log-likelihood
            kresids <- y - kcontrib_pop - kcontrib_sub
            samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

            ## calculate unnormalised log-posterior
            para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
            samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
                samples$ll[k]
        }
    }

    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}

## This is a Gibbs sampler v2 for longitudinal Bayesian ridge; see thesis.
## Update two block of parameters: variance and coefs
## Requirements: Bmat, y, grp, Kmat, dim_sub1
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_sub_v2 <- function(y, grp, Bmat, Kmat, dim_sub1, burn, size,
                               init = NULL, prior = NULL, verbose = TRUE) {

    ## hyperparemeters for priors
    if (is.null(prior)) {
        prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      iw_v = dim_sub1 + 1,
                      iw_lambda = diag(dim_sub1)) # assuming symmetry
    }

    check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))
    if (check_names == 0) {
        stop("Missing prior hyperparameters.")
    }

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    idx <- tapply(seq_len(n_samples), grp, function(x) x)

    ## some crossproducts precalculation
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    Bxy_sub <- matrix(NA, n_terms, n_subs, dimnames = list(NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
        Bxy_sub[, i] <- crossprod(Bmat[idx[[i]], ], y[idx[[i]]])
    }
    xKmat <- crossprod(Kmat)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)),
                    lp = rep(NA, size),
                    ll = rep(NA, size))

    ## initialise theta with a penalised LS estimate, and delta with rnorms with
    ## small sd
    pls <- tcrossprod(solve(crossprod(Bmat) + xKmat), Bmat) %*% as.vector(y)
    init <- initialise_with_pls(init, n_terms, grp, pls)
    kcoef <- list(pop = init$pop, sub = init$sub)

    ## initialise some intermediate output
    BMB <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    BMy <- matrix(NA, n_terms, n_subs, dimnames = list(NULL, levels(grp)))

    ## initialise prediction contribution by population coefs and subjects deviations
    kcontrib_pop <- Bmat %*% kcoef$pop
    kcontrib_sub <- rep(NA, n_samples)
    for (i in levels(grp)) {
        kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
    }

    for (k in seq.int(-burn + 1, size)) {
        ## update precisions
        para <- list(Kmat = Kmat, rank_K = rank_K, dim_sub1 = dim_sub1,
                     n_samples = n_samples, n_subs = n_subs, n_terms = n_terms)
        kresids <- y - kcontrib_pop - kcontrib_sub
        kprec <- update_prec(kcoef, kresids, para, prior)

        kprec_sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms - dim_sub1))

        ## update theta
        for (i in levels(grp)) {
            ## for numerical stability, these steps are simplified
            xBmat_i <- xBmat_sub[, , i]
            Li <- xBmat_i + kprec_sub / kprec$eps
            inv_Li <- chol2inv(chol(Li))
            BMB[, , i] <- kprec$eps * (diag(n_terms) - xBmat_i %*% inv_Li) %*% xBmat_i
            ## BMB[, , i] <- kprec$eps * xBmat_i - xBmat_i %*% inv_Li %*% xBmat_i

            Bxy_i <- Bxy_sub[, i]
            BMy[, i] <- kprec$eps * (Bxy_i - xBmat_i %*% inv_Li %*% Bxy_i)
            ## BMy[, i] <- kprec$eps * Bxy_i - xBmat_i %*% inv_Li %*% Bxy_i
        }
        Phi <- kprec$pop * xKmat + rowSums(BMB, dims = 2)
        inv_Phi <- chol2inv(chol(Phi))
        kcoef$pop <- t(mvtnorm::rmvnorm(1, inv_Phi %*% rowSums(BMy), inv_Phi))
        kcontrib_pop <- Bmat %*% kcoef$pop

        ## update delta
        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec$eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec$eps
            kcoef$sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec$pop
            samples$precision$sub1[, , k] <- kprec$sub1
            samples$precision$sub2[k] <- kprec$sub2
            samples$precision$eps[k] <- kprec$eps
            samples$population[, k] <- kcoef$pop
            samples$subjects[, , k] <- kcoef$sub

            ## calculate unnormalised log-likelihood
            kresids <- y - kcontrib_pop - kcontrib_sub
            samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

            ## calculate unnormalised log-posterior
            para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
            samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
                samples$ll[k]
        }
    }
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}

## initialise theta with constraint, and delta with tnorms with small sd
initialise_with_Amat <- function(init, n_terms, grp, Amat) {
    ## init: NULL or list(pop, sub); any of the element can be NULL

    Ainv <- diag(NCOL(Amat))
    Ainv[row(Ainv) > diff(dim(Amat))] <- Amat
    Ainv <- solve(Ainv)
    n_subs <- length(unique(grp))

    if (is.null(init$pop)) {
        init$pop <- Ainv %*% rep(1, ncol(Amat))
    } else {
        if (length(init$pop) == n_terms) {
            init$pop <- as.vector(init$pop)
            cat("Population initial values supplied.\n")
        } else {
            stop("Invalid dimension of population initial values.")
        }
    }
    if (is.null(init$sub)) {
        ## this is coming from gen_init function
        lower_right <- (-Amat %*% init$pop) + 1
        lower_left <- init$pop[1:diff(dim(Ainv))]
        init$sub <- t(tnorm::rmvtnorm(n_subs, init$pop, diag(n_terms) * 0.01,
                                       initial = Ainv %*% c(lower_left, lower_right),
                                       F = Amat,
                                       g = Amat %*% init$pop))
        colnames(init$sub) <- levels(grp)
    } else {
        if (dim(init$sub) == c(n_terms, n_subs)) {
            init$sub <- as.matrix(init$sub)
            cat("Subjects initial values supplied.\n")
            if (is.null(colnames(init$sub))) {
                colnames(init$sub) <- levels(grp)
            } else {
                init$sub <- init$sub[, levels(grp)]
            }
        } else {
            stop("Invalid dimension of subject initial values.")
        }
    }

    ## check feasibility
    if (any(Amat %*% init$pop < 0)) {
        stop("Population initial value violates constraints.")
    }
    if (any(Amat %*% (as.numeric(init$pop) + init$sub) < 0)) {
        stop("Subjects initial value violates constraints.")
    }

    init$Ainv <- Ainv
    init
}

## This is a Gibbs sampler for constrained longitudinal Bayesian ridge; see thesis.
## A * theta >= 0 ; A * (theta + delta_i) >= 0
## The vector of constriant has to be zero.
## Assumption: Full row rank for Amat.
## Requirements: Bmat, y, grp, Kmat, dim_sub1, Amat
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_cons_sub <- function(y, grp, Bmat, Kmat, dim_sub1, Amat, burn, size,
                                 init = NULL, prior = NULL, prec = NULL,
                                 verbose = TRUE) {

    ## hyperparemeters for priors
    if (is.null(prior)) {
        prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      iw_v = dim_sub1 + 1,
                      iw_lambda = diag(dim_sub1)) # assuming symmetry
    }

    check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))
    if (check_names == 0) {
        stop("Missing prior hyperparameters.")
    }

    if (NCOL(Amat) != NCOL(Bmat)) {
        stop("Dims of Amat and Bmat inconsistent.")
    }

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    idx <- tapply(seq_len(n_samples), grp, function(x) x)

    ## some crossproducts precalculation
    xBmat <- crossprod(Bmat)
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
    }
    xKmat <- crossprod(Kmat)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)),
                    lp = rep(NA, size),
                    ll = rep(NA, size))

    ## initialise theta and delta with tnorms with small sd, and a matrix to
    ## generate feasible starting points quickly
    init <- initialise_with_Amat(init, n_terms, grp, Amat)
    kcoef <- list(pop = init$pop, sub = init$sub)
    Ainv <- init$Ainv
    gen_init <- function(lower, mu) {
        Ainv %*% c(mu[1:(NCOL(Ainv) - length(lower))], lower + 1)
    }

    ## initialise prediction contribution by population coefs and subjects deviations
    kcontrib_pop <- Bmat %*% kcoef$pop
    kcontrib_sub <- rep(NA, n_samples)
    for (i in levels(grp)) {
        kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
    }

    for (k in seq.int(-burn + 1, size)) {
        ## update precisions
        para <- list(Kmat = Kmat, rank_K = rank_K, dim_sub1 = dim_sub1,
                     n_samples = n_samples, n_subs = n_subs, n_terms = n_terms)
        kresids <- y - kcontrib_pop - kcontrib_sub
        kprec <- update_prec(kcoef, kresids, para, prior)

        if (!is.null(prec)) {
            kprec$pop <- prec$pop
            kprec$sub1 <- prec$sub1
            kprec$sub2 <- prec$sub2
        }

        ## update theta
        lower_pop <- apply(-Amat %*% kcoef$sub, 1, max)
        lower_pop[lower_pop < 0] <- 0
        M_pop <- chol2inv(chol(xBmat + kprec$pop / kprec$eps * xKmat))
        mu <- M_pop %*% crossprod(Bmat, y - kcontrib_sub)
        sig <- 1 / kprec$eps * M_pop
        ## kcoef$pop <- t(mvtnorm::rmvnorm(1, mu, sig))
        kcoef$pop <- t(tnorm::rmvtnorm(1, mu, sig,
                                       initial = gen_init(lower_pop, mu),
                                       ## initial = kcoef$pop,
                                       F = Amat,
                                       g = -lower_pop))
        kcontrib_pop <- Bmat %*% kcoef$pop

        ## update delta
        lower_sub <- -Amat %*% kcoef$pop
        kprec_sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms - dim_sub1))

        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec$eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec$eps
            ## kcoef$sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
            kcoef$sub[, i] <- t(tnorm::rmvtnorm(1, mu, sig,
                                                initial = gen_init(lower_sub, mu),
                                                ## initial = kcoef$sub[, i],
                                                F = Amat,
                                                g = -lower_sub))
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec$pop
            samples$precision$sub1[, , k] <- kprec$sub1
            samples$precision$sub2[k] <- kprec$sub2
            samples$precision$eps[k] <- kprec$eps
            samples$population[, k] <- kcoef$pop
            samples$subjects[, , k] <- kcoef$sub

            ## calculate unnormalised log-likelihood
            kresids <- y - kcontrib_pop - kcontrib_sub
            samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

            ## calculate unnormalised log-posterior
            para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
            samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
                samples$ll[k]

        }
    }

    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}

## This is a Gibbs sampler for constrained longitudinal Bayesian ridge, with an
## HMC step to generate pop and sub coefs; see thesis.
## A * theta >= 0 ; A * (theta + delta_i) >= 0
## The vector of constriant has to be zero.
## Assumption: Full row rank for Amat.
## Requirements: Bmat, y, grp, Kmat, dim_sub1, Amat
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_cons_sub_v2 <- function(y, grp, Bmat, Kmat, dim_sub1, Amat, burn, size,
                                    init = NULL, prior = NULL, prec = NULL,
                                    verbose = TRUE) {

    ## hyperparemeters for priors
    if (is.null(prior)) {
        prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      iw_v = dim_sub1 + 1,
                      iw_lambda = diag(dim_sub1)) # assuming symmetry
    }

    check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))
    if (check_names == 0) {
        stop("Missing prior hyperparameters.")
    }

    if (NCOL(Amat) != NCOL(Bmat)) {
        stop("Dims of Amat and Bmat inconsistent.")
    }

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    idx <- tapply(seq_len(n_samples), grp, function(x) x, simplify = FALSE)

    ## construct full constraint matrix
    fAmat_left <- t(matrix(t(Amat), NCOL(Amat), NROW(Amat) * (n_subs + 1)))
    fAmat_topright <- matrix(0, NROW(Amat), n_terms * n_subs)
    fAmat_btmright <- block_diag(Amat, size = n_subs)
    fAmat <- cbind(fAmat_left, rbind(fAmat_topright, fAmat_btmright))
    lower <- rep(0, NROW(Amat))
    flower <- rep(lower, n_subs + 1)

    ## construct full design matrix
    fBmat_right <- do.call(block_diag, lapply(idx, function(x) Bmat[x, ]))
    fBmat <- cbind(Bmat, fBmat_right)

    ## some crossproducts precalculation
    xfBmat <- crossprod(fBmat)
    fBxy <- crossprod(fBmat, y)
    xKmat <- crossprod(Kmat)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)),
                    lp = rep(NA, size),
                    ll = rep(NA, size))

    ## initialise theta and delta with tnorms with small sd
    init <- initialise_with_Amat(init, n_terms, grp, Amat)
    kcoef <- list(pop = init$pop, sub = init$sub)
    fkcoef <- c(kcoef$pop, kcoef$sub)

    ## initialise prediction contribution by population coefs and subjects deviations
    kcontrib_pop <- Bmat %*% kcoef$pop
    kcontrib_sub <- rep(NA, n_samples)
    for (i in levels(grp)) {
        kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
    }

    for (k in seq.int(-burn + 1, size)) {
        ## update precisions
        para <- list(Kmat = Kmat, rank_K = rank_K, dim_sub1 = dim_sub1,
                     n_samples = n_samples, n_subs = n_subs, n_terms = n_terms)
        kresids <- y - kcontrib_pop - kcontrib_sub
        kprec <- update_prec(kcoef, kresids, para, prior)

        ## fix precision if provided
        if (!is.null(prec)) {
            kprec$pop <- prec$pop
            kprec$sub1 <- prec$sub1
            kprec$sub2 <- prec$sub2
        }

        ## construct Sigma
        kprec_sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms - dim_sub1))
        fkprec_all <- block_diag(kprec_sub, size = n_subs + 1)
        fkprec_all[1:n_terms, 1:n_terms] <- kprec$pop * xKmat

        ## update theta and delta
        M <- chol2inv(chol(xfBmat + fkprec_all / kprec$eps))
        mu <- M %*% fBxy
        sig <- 1 / kprec$eps * M
        fkcoef <- t(tnorm::rmvtnorm(1, mu, sig, fkcoef, fAmat, -flower, 10))
        ## fkcoef <- t(mvtnorm::rmvnorm(1, mu, sig))  # no problem in unconstrained case
        kcoef$pop <- fkcoef[1:n_terms]
        kcoef$sub <- matrix(fkcoef[-(1:n_terms)], n_terms, n_subs,
                            dimnames = list(NULL, levels(grp)))
        kcontrib_pop <- Bmat %*% kcoef$pop
        for (i in levels(grp)) {
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec$pop
            samples$precision$sub1[, , k] <- kprec$sub1
            samples$precision$sub2[k] <- kprec$sub2
            samples$precision$eps[k] <- kprec$eps
            samples$population[, k] <- kcoef$pop
            samples$subjects[, , k] <- kcoef$sub

            ## calculate unnormalised log-likelihood
            kresids <- y - kcontrib_pop - kcontrib_sub
            samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

            ## calculate unnormalised log-posterior
            para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
            samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
                samples$ll[k]
        }
    }

    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}

