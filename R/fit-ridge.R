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
  fm$data <- data[c('x', 'y', 'sub', 'pop')]
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
#' @export
fit_bs_splines <- function(data, K, deg, size, burn, ridge = FALSE, init = NULL,
                           prior = NULL) {

  check_data_K(data, K)
  if (!("pop" %in% names(data))) data$pop <- "dme__"


  ## design matrix for population curves
  des_info_pop <- get_design_bs(data$x, K$pop, deg)
  des_info_sub <- get_design_bs(data$x, K$sub, deg)
  if (ridge) {
    ## ridge regression
    Bmat <- list(pop = des_info_pop$design, sub = des_info_sub$design)
    Kmat <- get_diff_mat(K$pop + deg + 1, deg + 1) # difference matrix
  } else {
    ## LMM
    n_bsf_pop <- K$pop + deg + 1    # number of basis functions
    n_bsf_sub <- K$sub + deg + 1    # number of basis functions
    D_pop <- get_diff_mat(K$pop + deg + 1, deg + 1) # difference matrix
    D_sub <- get_diff_mat(K$sub + deg + 1, deg + 1) # difference matrix
    Tmat_pop <- cbind(-1/sqrt(n_bsf_pop), poly(1:n_bsf_pop, deg = deg, raw = FALSE),
                      crossprod(D_pop, solve(tcrossprod(D_pop))))
    Tmat_sub <- cbind(-1/sqrt(n_bsf_sub), poly(1:n_bsf_sub, deg = deg, raw = FALSE),
                      crossprod(D_sub, solve(tcrossprod(D_sub))))
    Bmat <- list(pop = des_info_pop$design %*% Tmat_pop,
                 sub = des_info_sub$design %*% Tmat_sub)
    Kmat <- cbind(matrix(0, K$pop, deg + 1), diag(K$pop))
  }

  fm <- bayes_ridge_sub_v2(y = data$y, grp = list(pop = data$pop, sub = data$sub),
                           Bmat = Bmat, Kmat = Kmat, dim_sub1 = deg + 1, burn = burn,
                           size = size, init = init, prior = prior)

  fm$basis <- list(pop = NA, sub = NA)
  if (ridge) {
    fm$basis$pop <- list(type = 'bs-ridge', knots = des_info_pop$knots, degree = deg)
    fm$basis$sub <- list(type = 'bs-ridge', knots = des_info_sub$knots, degree = deg)
  } else {
    fm$basis$pop <- list(type = 'bs', trans_mat = Tmat_pop, knots = des_info_pop$knots,
                         degree = deg)
    fm$basis$sub <- list(type = 'bs', trans_mat = Tmat_sub, knots = des_info_sub$knots,
                         degree = deg)
  }
  fm$data <- data[c('x', 'y', 'sub', 'pop')]
  fm$data$pop <- as.factor(fm$data$pop)
  fm$data$sub <- as.factor(fm$data$sub)
  ## fm$data <- dplyr::select(data, .data$x, .data$y, .data$sub, .data$pop) %>%
  ##   dplyr::mutate(sub = as.factor(.data$sub), pop = as.factor(.data$pop))
  fm
}

}

## predict with posterior means
## need a function to do posterior predictive check
## predict.fsso <- function(object, newdata = NULL) {
  
##   if (is.null(newdata)) {
##     Bmat <- object$Bmat
##   } else {
##     stopifnot(is.data.frame(newdata),
##               'x' %in% names(newdata))
##     Bmat_pop <- get_model_mat(newdata$x, object$basis$pop)
##     Bmat_sub <- get_model_mat(newdata$x, object$basis$sub)
    
##   }

## }




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
#' @param K A list of numbers of interior knots for the population and subject
#'     curves. The list must contain \code{sub} and \code{pop}.
#' @param deg The degree of the spline polynomial
#' @param Xmat A design matrix for the non-spline terms, i.e. fixed and random
#'   effects.
#' @param ranef A numeric vector of indices of the Xmat columns where the
#'   corresponding coefficients should be treated as random effects. NULL if
#'   everything is fixed. This option will get overridden by the 'beta' term in
#'   `prec`.
#' @param size The number of samples to be drawn from the posterior.
#' @param burn The number of samples to burn before recording. Default to
#'     one-tenth of \code{size}.
#' @param ridge Use ridge regression or linear mixed effect models?
#' @param init List of matrices of the coefficients, containing \code{sub} and
#'     \code{pop}. The columnns of the matrices correspond to each
#'     population/subject.
#' @param prior List of hyperparameters of the covariance priors.
#'
#' @return A list with posterior means, samples and information of the basis
#'     functions.
#'
#' @export
fit_bs_splines_v3 <- function(data, K, deg, Xmat, ranef, size, burn,
                              ridge = FALSE, init = NULL, prior = NULL) {

  check_data_K(data, K)
  if (!("pop" %in% names(data))) data$pop <- "dme__"


  ## design matrix for population curves
  des_info_pop <- get_design_bs(data$x, K$pop, deg)
  des_info_sub <- get_design_bs(data$x, K$sub, deg)
  if (ridge) {
    ## ridge regression
    Bmat <- list(pop = des_info_pop$design, sub = des_info_sub$design)
    Kmat <- get_diff_mat(K$pop + deg + 1, deg + 1) # difference matrix
  } else {
    ## LMM
    n_bsf_pop <- K$pop + deg + 1    # number of basis functions
    n_bsf_sub <- K$sub + deg + 1    # number of basis functions
    D_pop <- get_diff_mat(K$pop + deg + 1, deg + 1) # difference matrix
    D_sub <- get_diff_mat(K$sub + deg + 1, deg + 1) # difference matrix
    Tmat_pop <- cbind(-1/sqrt(n_bsf_pop), poly(1:n_bsf_pop, deg = deg, raw = FALSE),
                      crossprod(D_pop, solve(tcrossprod(D_pop))))
    Tmat_sub <- cbind(-1/sqrt(n_bsf_sub), poly(1:n_bsf_sub, deg = deg, raw = FALSE),
                      crossprod(D_sub, solve(tcrossprod(D_sub))))
    Bmat <- list(pop = des_info_pop$design %*% Tmat_pop,
                 sub = des_info_sub$design %*% Tmat_sub)
    Kmat <- cbind(matrix(0, K$pop, deg + 1), diag(K$pop))
  }

  if (is.null(prior)) {
    prior <- list(theta = list(a = 0.001, b = 0.001),
                  delta = list(a = 0.001, b = 0.001,
                               v = deg + 2, lambda = diag(deg + 1)),
                  ## beta = list(a = 0.001, b = 0.001),
                  eps = list(a = 0.001, b = 0.001))
  }
  
  fm <- bayes_ridge_semi(y = data$y, grp = list(pop = data$pop, sub = data$sub),
                         Bmat = Bmat, Xmat = Xmat, Kmat = Kmat, dim_block = deg + 1,
                         ranef = ranef, burn = burn, size = size, init = init,
                         prior = prior, prec = NULL)

  fm$basis <- list(pop = NA, sub = NA)
  if (ridge) {
    fm$basis$pop <- list(type = 'bs-ridge', knots = des_info_pop$knots, degree = deg)
    fm$basis$sub <- list(type = 'bs-ridge', knots = des_info_sub$knots, degree = deg)
  } else {
    fm$basis$pop <- list(type = 'bs', trans_mat = Tmat_pop, knots = des_info_pop$knots,
                         degree = deg)
    fm$basis$sub <- list(type = 'bs', trans_mat = Tmat_sub, knots = des_info_sub$knots,
                         degree = deg)
  }

  fm$data <- data[c('x', 'y', 'sub', 'pop')]
  fm$data$pop <- as.factor(fm$data$pop)
  fm$data$sub <- as.factor(fm$data$sub)

  fm$model_mat <- Bmat
  fm$pop_of_subs <- tapply(data$pop, data$sub,
                           function(x) as.character(unique(x)),
                           simplify = FALSE)

  ## fm$coefficients <- 
  ## fm$data <- dplyr::select(data, .data$x, .data$y, .data$sub, .data$pop) %>%
  ##   dplyr::mutate(sub = as.factor(.data$sub), pop = as.factor(.data$pop))
  class(fm) <- 'fsso'
  fm
}





