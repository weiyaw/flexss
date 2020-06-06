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
#' @param spline A list of formulae of splines s(x, by =, K = , deg = )
#' @param fixed A formula ~x1 + x2 ...
#' @param random A formula ~u1 + u2 ...
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
fit_bs_splines_v3 <- function(data, K, deg, spline, fixed = NULL, random = NULL,
                              Xmat = NULL, ranef = NULL,
                              size, burn, ridge = FALSE, init = NULL, prior = NULL) {

  check_data_K(data, K)
  if (!("pop" %in% names(data))) data$pop <- "dme__"

  ## make grp$sub to stick together. sorting the data according to the subject
  ## index is the quickest way.
  ## grp$sub should be nested within grp$pop.
  ## the grp$pop are also sticked together, because this is potentially useful
  fit_data <- data[order(data$pop, data$sub), ]
  
  ## design matrix for population curves
  dmatinfo_pop <- get_design_bs(fit_data$x, K$pop, deg)
  dmatinfo_sub <- get_design_bs(fit_data$x, K$sub, deg)
  if (ridge) {
    ## ridge regression
    Bmat <- list(pop = dmatinfo_pop$design, sub = dmatinfo_sub$design)
    Kmat <- get_diff_mat(K$pop + deg + 1, deg + 1) # difference matrix
  } else {
    ## LMM
    ## Tmat <- list(pop = get_transform_bs(K$pop + deg + 1, deg + 1),
    ##              sub = get_transform_bs(K$sub + deg + 1, deg + 1))
    ## Bmat <- list(pop = dmatinfo_pop$design %*% Tmat$pop,
    ##              sub = dmatinfo_sub$design %*% Tmat$sub)
    Kmat <- cbind(matrix(0, K$pop, deg + 1), diag(K$pop))
  }

  spline_mat <- purrr::map(spline, parse_spline, data = fit_data)
  fixed_mat <- parse_effect(fixed, fit_data)
  random_mat <- parse_effect(random, fit_data)

  ## this chunck eventually to be replaced when additive model is done
  Tmat <- list(pop = spline_mat$pop$trans_mat,
               sub = spline_mat$sub$trans_mat)
  Bmat <- list(pop = spline_mat$pop$model_mat,
               sub = spline_mat$sub$model_mat)
  Xmat <- fixed_mat$model_mat
  Kmat <- Kmat[, -1]
  
  if (is.null(prior)) {
    prior <- list(theta = list(a = 0.01, b = 0.01),
                  delta = list(a = 0.01, b = 0.01,
                               v = deg + 2, lambda = diag(deg + 1)),
                  ## beta = list(a = 0.001, b = 0.001),
                  eps = list(a = 0.01, b = 0.01))
  }

  grp <- purrr::map(fit_data[c('pop', 'sub')], ~factor(.x, levels = unique(.x)))
  ## grp_pop_names <- setNames(unique(grp$pop), unique(grp$pop))
  pop_names <- unique(grp$pop)
  Bmat_pop_ls <- purrr::map(pop_names, ~(Bmat$pop * (grp$pop == .x)))
  Bmat$pop_expand <- do.call(cbind, Bmat_pop_ls)
  colnames(Bmat$pop_expand) <- rep(pop_names, each = NCOL(Bmat$pop))
  Kmat_expand <- block_diag(Kmat, size = length(Bmat_pop_ls))  
  
  ## check if cbind(Bmat_pop, Xmat) is full column rank
  stopifnot(det(crossprod(cbind(Bmat$pop_expand, Xmat))) > 0)
  
  fm <- bayes_ridge_semi(y = fit_data$y,
                         grp = grp,
                         Bmat = list(pop = Bmat$pop_expand, sub = Bmat$sub),
                         Xmat = Xmat,
                         Kmat = Kmat_expand,
                         dim_block = deg + 1,
                         ranef = ranef, burn = burn, size = size, init = init,
                         prior = prior, prec = NULL)

  
  ## convert pop coef back to array 
  dim(fm$samples$coef$theta) <- c(NCOL(Bmat$pop), length(pop_names), size)
  dimnames(fm$samples$coef$theta) <- list(NULL, pop_names, NULL)
  dim(fm$means$coef$theta) <- c(NCOL(Bmat$pop), length(pop_names))
  dimnames(fm$means$coef$theta) <- list(NULL, pop_names)

  
  fm$basis <- list(pop = NA, sub = NA)
  if (ridge) {
    fm$basis$pop <- list(type = 'bs-ridge',
                         knots = dmatinfo_pop$knots,
                         degree = deg)
    fm$basis$sub <- list(type = 'bs-ridge',
                         knots = dmatinfo_sub$knots,
                         degree = deg)
  } else {
    fm$basis$pop <- list(type = 'bs', trans_mat = Tmat$pop,
                         knots = dmatinfo_pop$knots,
                         degree = deg)
    fm$basis$sub <- list(type = 'bs', trans_mat = Tmat$sub,
                         knots = dmatinfo_sub$knots,
                         degree = deg)
  }

  ## return the fit_data, not the original dataset
  fm$data <- data
  ## fm$data$pop <- as.factor(fm$data$pop)
  ## fm$data$sub <- as.factor(fm$data$sub)

  ## fm$spline <- spline_mat
  ## fm$fixed <- fixed_mat
  ## fm$random <- random_mat
  fm$spline <- purrr::map(spline_mat, ~`[<-`(.x, 'model_mat', NULL))
  fm$fixed <- `[<-`(fixed_mat, 'model_mat', NULL)
  fm$random <- `[<-`(random_mat, 'model_mat', NULL)
  ## fm$formula <- list(spline = spline,
  ##                    fixed = fixed,
  ##                    random = random)

  ## fm$model_mat <- list(Bmat = Bmat, Xmat = Xmat)
  fm$pop_of_subs <- tapply(fit_data$pop, fit_data$sub,
                           function(x) as.character(unique(x)),
                           simplify = FALSE)

  ## fm$coefficients <- 
  ## fm$data <- dplyr::select(data, .data$x, .data$y, .data$sub, .data$pop) %>%
  ##   dplyr::mutate(sub = as.factor(.data$sub), pop = as.factor(.data$pop))
  class(fm) <- 'fsso'
  fm
}


## ff <- y ~ s(x, by = pop) + s(x, by = sub)


## transform ~s(x, by = , K = , deg =) to bs model mat (and associated details)
## everything is evaluted at the scope where the formula is defined
## return: list(Bmat, Tmat, by, K, knotsdeg)
parse_spline <- function(fo, data) {
  if (class(fo) == 'formula' && length(fo) == 2) {
    qo <- fo[[2]] # qo = quote
  } else {
    stop("spline not specified as a one-sided formula.")
  }

  if (qo[[1]] != 's') {
    stop("formula is not 's'.")
  }
  
  res <- list(formula = fo)
  ## response and by variables are evaluated on data
  x <- eval(qo[[2]], data)
  res$by <- eval(qo[['by']], data)

  ## the rest of argument evaluated on the env where it was defined
  stopifnot(!is.null(qo[['K']]), !is.null(qo[['deg']]))
  K <- eval(qo[['K']], attr(fo, '.Environment'))
  intercept <- eval(qo[['intercept']], attr(fo, '.Environment'))
  res$degree <- eval(qo[['deg']], attr(fo, '.Environment'))

  dmatinfo <- get_design_bs(x, K, res$degree)
  if (!is.null(intercept) && intercept) {
    res$trans_mat <- get_transform_bs(K + res$degree + 1, res$degree + 1)
  } else {
    ## remove the intercept of the spline
    res$trans_mat <- get_transform_bs(K + res$degree + 1, res$degree + 1)[, -1]
  }
  res$model_mat <- dmatinfo$design %*% res$trans_mat
  res$type <- 'bs'
  res$knots <- dmatinfo$knots
  res
}

parse_effect <- function(fo, data) {
  if (is.null(fo)) {
    NULL
  } else {
    ## mt: model term
    mt <- terms(fo, data = data)
    stopifnot(attr(mt, 'response') == 0) # must not have response term in the formula

    variable <- attr(terms(fo), 'variable')
    attr(mt, 'dataClasses') <- eval(variable, purrr::map(data, .MFclass))
    mode(attr(mt, 'dataClasses')) <- 'character'
    names(attr(mt, 'dataClasses')) <- as.character(variable[-1])

    ## mm: model matrix
    mm <- model.matrix(mt, data)
    list(terms = mt,
         model_mat = mm,
         assign = attr(mm, 'assign'),
         contrasts = attr(mm, 'contrasts'),
         xlevels = .getXlevels(mt, data))
    ## attr(res$terms, 'dataClasses')
  }
}



## efficient prediction calculation of factor variables
## model_mat: the model matrix to be multiplied
## grp: a vector specifying the group which each row of model matrix corresponds to
## coef: a matrix of spline coef. colnames must have unique(grp)
predict_grp <- function(model_mat, grp, coef) {

  stopifnot(unique(grp) %in% colnames(coef),
            length(grp) == NROW(model_mat))

  res <- rep(NA, length.out = length(grp))

  for (l in unique(grp)) {
    idx <- grp == l
    res[idx] <- model_mat[idx, ] %*% coef[, l]
  }
  res
}




## predict with posterior means
## need a function to do posterior predictive check
predict.fsso <- function(object, newdata = NULL) {

  ## temporary measure. change the output samples in the future
  names(object$means$coef) <- c('pop', 'sub', 'beta')
    
  if (is.null(newdata)) {
    newdata <- object$data
  }  

  predict_spl_i <- function(spl, spl_name) {
    x <- eval(spl$formula[[2]][[2]], newdata)
    model_mat <- get_model_mat(x, spl)
    grp <- eval(spl$formula[[2]][['by']], newdata)
    predict_grp(model_mat, grp, object$means$coef[[spl_name]])
  }
  
  res1 <- purrr::imap(object$spline, predict_spl_i)

  predict_eff_i <- function(effect, newdata, beta) {
    if (is.null(effect)) {
      0
    } else {
      mt <- terms(effect)
      mf <- model.frame(mt, newdata, xlev = effect$xlevels)
      .checkMFClasses(attr(mt, 'dataClasses'), mf)
      mm <- model.matrix(mt, mf, contrast.arg = effect$contrasts)
      stopifnot(colnames(mm) == names(beta))
      as.numeric(mm %*% beta)
    }
  }

  res2 <- list(fixed = predict_effect(object[['fixed']], newdata, object$means$coef$beta))

  c(res1, res2)
  ## purrr::reduce(c(res1, res2), `+`)

}


