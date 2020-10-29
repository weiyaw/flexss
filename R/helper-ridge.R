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
            message("Population initial values supplied.")
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
            message("Subjects initial values supplied.")
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


#' Construct a difference matrix
#'
#' \code{get_diff_mat} returns a difference matrix of an arbitrary order and size.
#'
#' This function returns a \eqn{k^{th}} order difference matrix \eqn{D}, such
#' that \eqn{Dx} gives a vector of \eqn{k^{th}} order differences of \eqn{x}.
#' The parameter \code{size} is usually the dimension of \eqn{x}.
#'
#' @param size the column size of the difference matrix. Required.
#'
#' @param k the order of difference. Required.
#'
#' @return A \code{size - k} by \code{size} matrix.

get_diff_mat <- function(size, k) {
    if (size <= k) {
        stop("Order of difference greater than column size.")
    }

    D <- diag(1, size)
    D[(col(D) + 1) == row(D)] <- -1
    D <- D[-1, ]
    if (k == 1) {
      D
    } else {
      Recall(size - 1, k - 1) %*% D
    }
}

## A transformation matrix to enforce differnce penalty of B-spline
## 
## Coefs of B-spline are usually penalised with a k-order difference
## penalty. This function returns a matrix that transform the coefs such that
## the difference penalty on the coefs can be equivalently enforced by
## penalising the first 'k' elements of the transformed coefs to 0. That is, if
## \theta is the B-spline coef and its kth-order difference is penalised, then
## it is equivalent to shrink the first 'k' elements of \theta' to zero, where
## \theta = T %*% \theta'. This function returns the transformation matrix T.
##
## size: the number of coefficients/basis functions
## k: order of difference penalty
## return: a transformation matrix
get_transform_bs <- function(size, k) {
  Dmat <- get_diff_mat(size, k)
  Tmat <- cbind(-1/sqrt(size),
                stats::poly(seq_len(size), degree = k - 1, simple = TRUE),
                crossprod(Dmat, chol2inv(chol(tcrossprod(Dmat)))))
  unname(Tmat)
}

#' Construct an equidistant B-splines design matrix
#'
#' \code{get_design_bs} returns an equidistant B-splines design matrix without
#' explicitly specifying the location of the knots.
#'
#' The knots locations are the minimum and maximum of \code{x}, and \code{K}
#' equidistant points between these two extrema. The B-splines are equidistant,
#' i.e. no multiple knots at both ends of \eqn{x}. This function uses
#' \code{splines::splineDesign}.
#' 
#' @param x predictor vector. Required.
#' @param K number of inner knots, or a vector of all knots (interior,
#'     boundaries, and knots outside extrema). Required.
#' @param deg the degree of polynomial which the B-splines span. Required.
#' @param EPS tolerance error.
#'
#' @return a list with components `design' (\code{length(x)} by \code{K + deg +
#'   1} design matrix) and all `knots' (interior and boundaries knots, and knots
#'   outside extrema). If type == 'LMM', the transformation matrix is returned
#'   in the 'transform' field.
get_design_bs <- function(x, K, deg, EPS = 1e-6) {
    res <- list()

    ## get the knots
    if (length(K) == 1 && is.numeric(K)) {
        dist <- diff(range(x)) / (K + 1)
        knots <- seq(min(x) - (deg * dist) - EPS, max(x) + (deg * dist) + EPS,
                     len = K + 2 * (deg + 1))
        names(knots) <- NULL
        res$knots <- knots
    } else if (is.vector(K) && is.numeric(K)) {
        knots <- K
        res$knots <- knots
    } else {
        stop("Supplied knots must a numeric vector.")
    }
    names(knots) <- NULL

    res$design <- splines::splineDesign(knots, x, ord = deg + 1)
    res$knots <- knots
    res
}


## Return a design matrix at quartile (or supplied) knots, and all knots.
## The knots are taken at the quartiles, and min and max of x.
## x: all the explanatory data
## K: number of inner knots, or a vector of all knots (including extrema)
## deg: degree of the spline polynomial
get_design_tpf <- function(x, K, deg) {

    res <- list()

    ## get the inner knots
    if (length(K) == 1 && is.numeric(K)) {
        res$knots <- unname(stats::quantile(unique(x), seq(0, 1, len = K + 2)))
        ## res$knots <- unname(seq(min(x), max(x), length.out = K + 2))
        knots <- res$knots[-c(1, K + 2)]
    } else if (is.vector(K) && is.numeric(K)) {
        knots <- K[-c(1, length(K))]
        res$knots <- K
    } else {
        stop("Supplied knots must a numeric vector.")
    }
    names(knots) <- NULL

    ## get the design matrix
    splines <- outer(x, knots, `-`)

    ## this chunck can be generalised to higher degree polynomials
    if (deg > 0 && deg < 4) {
        splines <- splines^deg * (splines > 0)
        res$design <- cbind(1, stats::poly(x, degree = deg, raw = TRUE, simple = TRUE),
                            splines, deparse.level = 0)
        colnames(res$design) <- NULL
    } else {
        stop("Invalid degree.")
    }
    res
}


## Get a penalised least squares estimate with a given response, design matrix
## and penalty matrix.
get_pls <- function(response, design, penalty) {
  inv_term <- chol2inv(chol(crossprod(design) + crossprod(penalty)))
  tcrossprod(inv_term, design) %*% as.vector(response)
}

## get the maximum-a-posteriori or maximum likelihood estimate from the model
get_max <- function(fm, type = "map") {
    if (type == "map") {
        if (is.null(fm$samples$lp)) {
            stop("Log-posterior not available.")
        } else {
            idx <- which.max(fm$samples$lp)
        }
    } else if (type == "mle") {
        if (is.null(fm$samples$ll)) {
            stop("Log-likelihood not available.")
        } else {
            idx <- which.max(fm$samples$ll)
        }
    } else {
        stop("Unknown type in get_max.")
    }
    list(population = fm$samples$population[, idx],
         subjects = fm$samples$subjects[, , idx],
         precision = list(pop = fm$samples$precision$pop[idx],
                          sub1 = fm$samples$precision$sub1[, , idx],
                          sub2 = fm$samples$precision$sub2[idx],
                          eps = fm$samples$precision$eps[idx]))
}







