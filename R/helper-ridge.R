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

## initialise theta with a penalised LS estimate, and delta with rnorms with small sd
initialise_with_pls <- function(init, n_terms, grp, pls) {
    ## init: NULL or list(pop, sub); any of the element can be NULL

    n_pops <- length(unique(grp$pop))
    n_subs <- length(unique(grp$sub))
    if (is.null(init$pop)) {
        init$pop <- matrix(pls, n_terms, n_pops)
    } else {
        if (dim(init$sub) == c(n_terms, n_pops)) {
            init$pop <- as.matrix(init$sub)
            cat("Population initial values supplied.\n")
        } else {
            stop("Invalid dimension of population initial values.")
        }
    }
    if (is.null(init$sub)) {
        init$sub <- matrix(rnorm(n_terms * n_subs) * 0.01, n_terms, n_subs,
                           dimnames = list(NULL, levels(grp$sub)))
    } else {
        if (dim(init$sub) == c(n_terms, n_subs)) {
            init$sub <- as.matrix(init$sub)
            cat("Subjects initial values supplied.\n")
            if (is.null(colnames(init$sub))) {
                colnames(init$sub) <- levels(grp$sub)
            } else {
                init$sub <- init$sub[, levels(grp$sub)]
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

check_grp <- function(grp) {

    ## check grp field names
    if (all(c("pop", "sub") %in% names(grp))) {
        stop("Missing grp variables.")
    }

    ## check grp lengths
    if (length(grp$pop) != length(grp$sub)) {
        stop("Unequal length of grp$pop and grp$sub")
    }

    ## sub factor must be nested within pop factor
    sublvl_in_pop <- tapply(grp$sub, grp$pop, unique, simplify = FALSE)

    if (length(unlist(sublvl_in_pop)) > length(unique(grp$sub))) {
        stop("Each subject can only appear in one population.")
    }

    if (all(lapply(grp, is.factor))) {
        lapply(grp, droplevels)
    } else {
        lapply(grp, factor)
    }
}

check_Bmat <- function(Bmat, Kmat) {

    ## check grp field names
    if (all(c("pop", "sub") %in% names(Bmat))) {
        stop("Missing Bmat variables.")
    }

    ## check Bmat lengths
    if (NROW(Bmat$pop) != NROW(Bmat$sub)) {
        stop("Unequal length of Bmat$pop and Bmat$sub.")
    }

    ## check col of Bmat
    if (NCOL(Bmat$pop) != NCOL(Kmat)) {
        stop("Inconsistant columns of Bmat and Kmat.")
    }
    Bmat
}

check_prior <- function(prior, dim_sub1) {

    ## build a prior for covariance if not give
    if (is.null(prior)) {
        prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                      iw_v = dim_sub1 + 1,
                      iw_lambda = diag(dim_sub1)) # assuming symmetry
    }

    ## check prior field name
    check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                        c("pop", "sub2", "eps") %in% names(prior$ig_a),
                        c("pop", "sub2", "eps") %in% names(prior$ig_b))
    if (check_names == 0) {
        stop("Missing prior hyperparameters.")
    }
    prior
}

calc_xB <- function(Bmat, idx_sub, para) {

    ## crossproduct of Bmats (pop x pop, sub x pop, sub x sub)
    n_subs <- para$n_subs
    n_terms_pop <- para$n_terms_pop
    n_terms_sub <- para$n_terms_sub

    xB_pop <- array(NA, c(n_terms_pop, n_terms_pop, n_subs),
                    list(NULL, NULL, names(idx_sub)))
    xB_subpop <- array(NA, c(n_terms_sub, n_terms_pop, n_subs),
                       list(NULL, NULL, names(idx_sub)))
    xB_sub <- array(NA, c(n_terms_sub, n_terms_sub, n_subs),
                    list(NULL, NULL, names(idx_sub)))

    for (i in names(idx_sub)) {
        xB_pop[, , i] <- crossprod(Bmat$sub[idx_sub[[i]], ])
        xB_subpop[, , i] <- crossprod(Bmat$sub[idx_sub[[i]], ], Bmat$pop[idx_sub[[i]], ])
        xB_sub[, , i] <- crossprod(Bmat$pop[idx_sub[[i]], ])
    }
    list(pop = xB_pop, subpop = xB_subpop, sub = xB_sub)
}

calc_Bxy <- function(y, Bmat, idx_sub, para) {

    n_terms_sub <- para$n_terms_sub
    n_subs <- para$n_subs

    ## crossproduct of Bmat and y (Bmat_pop x y, Bmat_sub x y)
    Bxy_pop <- matrix(NA, n_terms_sub, n_subs, dimnames = list(NULL, names(idx_sub)))
    Bxy_sub <- matrix(NA, n_terms_sub, n_subs, dimnames = list(NULL, names(idx_sub)))

    for (i in names(idx_sub)) {
        Bxy_pop[, i] <- crossprod(Bmat$pop[idx_sub[[i]], ], y[idx_sub[[i]]])
        Bxy_sub[, i] <- crossprod(Bmat$sub[idx_sub[[i]], ], y[idx_sub[[i]]])
    }
    list(pop = Bxy_pop, sub = Bxy_sub)
}









