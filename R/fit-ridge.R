## check the validity of the dataset and K (knots)
check_data_K <- function(data, K) {
    if (all(c("y", "x", "sub") %in% names(data))) {
        if ("pop" %in% names(data)) {
            message("Fitting mean curves according to the populations.")
        } else {
            message("Fitting one mean curve for the whole dataset.")
        }
    } else {
        stop("'x', 'y', and 'sub' must exist in data.")
    }

    if (!all(c("sub", "pop") %in% names(K))) {
        stop("'K' must contain 'sub' and 'pop'.")
    }
}


#' @export
fit_tpf_splines <- function(data, K, deg, size, burn, init = NULL, prior = NULL) {

    check_data_K(data, K)
    if (!("pop" %in% names(data))) data$pop <- "dme__"

    ## design matrix for population curves
    des_info_pop <- get_design_tpf(data$x, K$pop, deg)
    Bmat <- list(pop = des_info_pop$design)
    Kmat <- cbind(matrix(0, K$pop, deg + 1), diag(K$pop))

    ## design matrix for subject-specific curves
    des_info_sub <- get_design_bs(data$x, K$sub, deg)
    Bmat$sub <- des_info_sub$design

    fm <- bayes_ridge_sub_v2(y = data$y, grp = list(pop = data$pop, sub = data$sub),
                             Bmat = Bmat, Kmat = Kmat, dim_sub1 = deg + 1, burn = burn,
                             size = size, init = init, prior = prior)

    fm$basis <- list(pop = NA, sub = NA)
    fm$basis$pop <- list(type = 'tpf', knots = des_info_pop$knots, degree = deg)
    fm$basis$sub <- list(type = 'tpf', knots = des_info_sub$knots, degree = deg)
    fm$data <- dplyr::select(data, x, y, sub, pop)
    fm
}

#' Fit a B-spline model with population and subject-specific curves.
#'
#' Fit a B-spline model which gives population and subject-specific curves.
#'
#' This model is fitted in a Bayesian framework, that is, this routine gives the
#' samples drawing from the posterior distribution associated to each regression
#' coefficients. The population (mean) curves are penalised with a difference
#' matrix of order \code{deg + 1}. The knots are equally spaced between
#' \code{range(data$x)}.
#'
#' @param data A data frame with at least three columns (\code{x}, \code{y} and
#'     \code{sub}). If a \code{pop} column is also given, separate mean curves
#'     are fitted to each population.
#'
#' @param K A list of numbers of interior knots for the population and subject
#'     curves. The list must contain \code{sub} and \code{pop}.
#'
#' @param deg The degree of the spline polynomial
#'
#' @param size The number of samples to be drawn from the posterior.
#'
#' @param burn The number of samples to burn before recording. Default to
#'     one-tenth of \code{size}.
#'
#' @param ridge Use ridge regression or linear mixed effect models?
#'
#' @param init List of matrices of the coefficients, containing \code{sub} and
#'     \code{pop}. The columnns of the matrices correspond to each
#'     population/subject.
#'
#' @param prior List of hyperparameters of the covariance priors.
#'
#' @return A list with posterior means, samples and information of the basis
#'     functions.
#'
fit_bs_splines <- function(data, K, deg, size, burn, ridge = FALSE, init = NULL,
                           prior = NULL) {

    check_data_K(data, K)
    if (!("pop" %in% names(data))) data$pop <- "dme__"

    n_bsf_pop <- K$pop + deg + 1    # number of basis functions
    D <- get_diff_mat(n_bsf_pop, deg + 1) # difference matrix

    ## design matrix for population curves
    des_info_pop <- get_design_bs(data$x, K$pop, deg)
    if (ridge) {
        ## ridge regression
        Bmat <- list(pop = des_info_pop$design)
        Kmat <- D
    } else {
        ## LMM
        Tmat <- cbind(-1/sqrt(n_bsf_pop), poly(1:n_bsf_pop, deg = deg, raw = FALSE),
                      crossprod(D, solve(tcrossprod(D))))
        Bmat <- list(pop = des_info_pop$design %*% Tmat)
        Kmat <- cbind(matrix(0, K$pop, deg + 1), diag(K$pop))
    }

    ## design matrix for subject-specific curves
    des_info_sub <- get_design_bs(data$x, K$sub, deg)
    Bmat$sub <- des_info_sub$design

    fm <- bayes_ridge_sub_v2(y = data$y, grp = list(pop = data$pop, sub = data$sub),
                             Bmat = Bmat, Kmat = Kmat, dim_sub1 = deg + 1, burn = burn,
                             size = size, init = init, prior = prior)

    fm$basis <- list(pop = NA, sub = NA)
    if (ridge) {
        fm$basis$pop <- list(type = 'bs-ridge', knots = des_info_pop$knots, degree = deg)
        fm$basis$sub <- list(type = 'bs-ridge', knots = des_info_sub$knots, degree = deg)
    } else {
        fm$basis$pop <- list(type = 'bs', trans_mat = Tmat, knots = des_info_pop$knots,
                             degree = deg)
        fm$basis$sub <- list(type = 'bs', trans_mat = Tmat, knots = des_info_sub$knots,
                             degree = deg)
    }
    fm$data <- dplyr::select(data, x, y, sub, pop)
    fm
}




