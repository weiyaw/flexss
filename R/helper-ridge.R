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


initialise_samples <- function(para, size, dim_sub1, grp) {
    n_terms_pop <- para$n_terms_pop
    n_terms_sub <- para$n_terms_sub
    n_subs <- para$n_subs
    n_pops <- para$n_pops
    list(population = array(NA, c(n_terms_pop, n_pops, size),
                            dimnames = list(NULL, levels(grp$pop))),
         subjects = array(NA, c(n_terms_sub, n_subs, size),
                          dimnames = list(NULL, levels(grp$sub))),
         precision = list(pop = rep(NA, size),
                          sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                          sub2 = rep(NA, size),
                          eps = rep(NA, size)),
         lp = rep(NA, size),
         ll = rep(NA, size))
}









