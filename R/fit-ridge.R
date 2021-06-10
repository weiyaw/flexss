## TO BE ABSORBED BY s(.)
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
#' The intercepts for all splines, other than the subject-specific splines, are
#' removed by default. This is to ensure that the intercept term in the spline
#' is not spanning the same space as the (potential) main effect terms (and
#' consequently having an unidentifiable model). For subject-specific splines,
#' all the coefficients are penalised, so including an intercept will not result
#' in an unidentifiable model.
#'
#'  Therefore, the intercept should be added manually either by setting
#' 'intercept = TRUE' or adding it as fixed effects. 
#' 
#' @param fixed A formula y ~ x1 + x2 ... specifying the response and fixed
#'   effects. This routine looks up the variable from the data, before looking
#'   at the environment where the formula is defined.
#' @param data A data frame with at least three columns (\code{x}, \code{y} and
#'   \code{sub}). If a \code{pop} column is also given, separate mean curves are
#'   fitted to each population.
#' @param spline A list of formulae of splines s(x, by =, ...). This routine
#'   looks up the variable from the data, before looking at the environment
#'   where the formula is defined. See [s()] for more details.
#' @param random A formula ~ u1 specifying the random effect term. Only one
#'   random effect is allowed at the moment. In the future, this argument might
#'   take a list of fomulae (so each formula corresponds to a random effect), or
#'   implement something like nlme::lme or lme4::lmer.
#' @param size The number of samples to be drawn from the posterior.
#' @param burn The number of samples to burn before recording. Default to
#'   one-tenth of \code{size}.
#' @param ridge Use ridge regression or linear mixed effect models?
#' @param init List of matrices of the coefficients, containing \code{sub} and
#'   \code{pop}. The columnns of the matrices correspond to each
#'   population/subject.
#' @param prior List of hyperparameters of the covariance priors.
#'
#' @return A list with posterior means, samples and information of the basis
#'     functions.
#'
#' @seealso [s()] for spline specification.
#' @export
fit_bs_splines_v4 <- function(fixed, data, spline, random = NULL,
                              size = 1000, burn = 0,
                              ridge = FALSE, init = NULL,
                              prior = NULL, prec = NULL) {

  ## check_data_K(data, K)
  ## if (!("pop" %in% names(data))) data$pop <- "dme__"

  ## make grp$sub to stick together. sorting the data according to the subject
  ## index is the quickest way.
  ## grp$sub should be nested within grp$pop.
  ## the grp$pop are also sticked together, because this is potentially useful
  spline <- purrr::map(spline, ~`[[<-`(.x, 2, match.call(s, .x[[2]])))
  which_sub <- which(purrr::map_lgl(spline, detect_sub))
  stopifnot(length(which_sub) == 1)
  fit_data <- data[reorder_data(spline[[which_sub]], data), ]
  
  spline_obj <- purrr::map(purrr::compact(spline), parse_spline, data = fit_data)
  effect_obj <- purrr::map(purrr::compact(list(fixed = fixed, random = random)),
                           parse_effect, data = fit_data)
  response <- parse_response(fixed, fit_data)
  ## LOOK AT FIXED EFFECT, NEED TO GIVE IT A NAME
  if (length(setdiff(names(spline_obj), '')) != length(spline_obj)) {
    names(spline_obj) <- paste0('spl', 1:length(spline_obj))
    message('Spline names overriden.')
  }

  fit_Bmat <- purrr::map(spline_obj, 'model_mat')
  fit_Xmat <- purrr::map(effect_obj, 'model_mat')
  
  ## check if fixed effects matrix is full column rank
  stopifnot(det(crossprod(fit_Xmat$fixed)) > 0)

  ## if no prec or prior is given
  if (is.null(prec) && is.null(prior)) {
    ## prior <- get_prior(fit_Bmat, fit_Xmat, a = 0.001, b = 0.001, v = 3)
    prior <- get_prior(fit_Bmat, fit_Xmat, a = -0.5, b = 0, v = -1)
  } else {
    ## prior <- get_prior(fit_Bmat, fit_Xmat, a = 0.001, b = 0.001, v = 3)
    prior <- get_prior(fit_Bmat, fit_Xmat, a = -0.5, b = 0, v = -1)
    message('For the time being, prior cannot be manually specified.')
  }
  
  fm <- bayes_ridge_semi_v4(y = response$y,
                            Bmat = fit_Bmat,
                            Xmat = fit_Xmat,
                            prior = prior,
                            prec = prec,
                            init = init,
                            burn = burn,
                            size = size)
  
  ## return the fit_data, not the original dataset
  fm$data <- data

  ## remove model_mat and return
  fm$spline <- purrr::map(spline_obj, ~`[<-`(.x, 'model_mat', NULL))
  fm$effect <- purrr::map(effect_obj, ~`[<-`(.x, 'model_mat', NULL))
  class(fm) <- 'blm' # bayesian longitudinal model (black lives matter)
  fm
}

## return an appropriate hyperparameter for the variance prior
get_prior <- function(Bmat, Xmat, a = -0.5, b = 0, v = -1, lambda = NULL) {
  res <- list()
  spl <- function(x) {
    if (attr(x, 'is_sub')) {
      lambda <- diag(attr(x, 'block_dim'))
      list(a = a, b = b, v = v, lambda = lambda)
    } else {
      list(a = a, b = b)
    }
  }

  ## spline hyperparameter
  res$spl <- purrr::map(Bmat, spl)
  eff <- function(x) {
    if (x == 'random') {
      list(a = a, b = b)
    } else {
      NULL
    }
  }
  
  ## effect hyperparameter
  if (is.null(Xmat)) {
    res$effect <- NULL
  } else {
    res$eff <- purrr::imap(Xmat, ~eff(.y))
  }
  res$eps <- list(a = a, b = b)
  res
}

## efficient prediction calculation of factor variables
## model_mat: the model matrix to be multiplied
## by: a vector specifying the group which each row of model matrix corresponds to
## coef: a matrix of spline coef. colnames must have unique(by)
predict_grp <- function(model_mat, by, coef) {
  
  if (is.null(by)) {
    stopifnot(NCOL(coef) == 1)
    as.numeric(model_mat %*% coef)
  } else {
    stopifnot(unique(by) %in% colnames(coef),
              length(by) == NROW(model_mat))
    res <- rep(NA, length.out = length(by))
    for (l in unique(by)) {
      idx <- by == l
      res[idx] <- model_mat[idx, ] %*% coef[, as.character(l)]
    }
    res
  }
}

## efficient prediction calculation of factor variables in bulk
## model_mat: the model matrix to be multiplied
## by: a vector specifying the group which each row of model matrix corresponds to
## coef: a 2D or 3D array of spline coef. The last dim are samples of coefs. If 3D, colnames must have unique(by).
## value: a 2D matrix, nrow = nrow of model_mat, ncol = number of samples
predict_grp2 <- function(model_mat, by, coefs) {
  
  if (is.null(by)) {
    stopifnot(length(dim(coefs)) == 2)
    model_mat %*% coefs
  } else {
    stopifnot(unique(by) %in% colnames(coefs),
              length(by) == NROW(model_mat))
    res <- matrix(NA, nrow = length(by), ncol = dim(coefs)[3])
    for (l in unique(by)) {
      idx <- by == l
      res[idx, ] <- model_mat[idx, ] %*% coefs[, as.character(l), ]
    }
    res
  }
}



#' Predict with posterior mean (or with other statistics)
#'
#' Return prediction with the posterior mean (or other statistics) of the
#' regression coefficients.
#' 
#' @param object an 'blm' object.
#' @param newdata a data frame with all the covariate required for prediction.
#' @param level a string specify the name of the spline. Only the component
#'   correponding to that spline is returned together with the effect terms).
#' @param fun a function applying to the samples of each marginal posterior of
#'   the regression coefs. By default, it's `mean`.
#' 
#' @export
predict.blm <- function(object, newdata = NULL, level = NULL, fun = NULL) {
  ## temporary measure. change the output samples in the future
  ## names(object$means$coef) <- c('pop', 'sub', 'fixed')
    
  if (is.null(newdata)) {
    newdata <- object$data
  }  

  predict_spl_i <- function(spl, coefs) {
    ## 'x' and 'by' variable must be in newdata, otherwise throw an error
    x <- eval(spl$call[['x']], newdata, baseenv())
    ## x <- newdata[[as.character(spl$call[['x']])]]
    model_mat <- get_model_mat(x, spl)
    grp <- eval(spl$call[['by']], newdata, baseenv())
    ## grp <- newdata[[as.character(spl$call[['by']])]]
    predict_grp(model_mat, grp, coefs)
  }

  ## make sure the levels are valid
  stopifnot(!is.null(names(object$spline)))
  if (is.null(level)) {
    level <- names(object$spline)
  }
  stopifnot(level %in% names(object$spline))

  ## predict with means if fun is NULL, or otherwise
  if (is.null(fun)) {
    coefs_ls <- object$means$coef
  } else {
    coefs_ls <- pstats_v4(object$samples$coef, fun = fun)
  }

  res1 <- purrr::map2(object$spline[level],
                      coefs_ls$spline[level],
                      predict_spl_i)


  predict_eff_i <- function(eff, coefs) {
    if (is.null(eff)) {
      0
    } else {
      mt <- stats::delete.response(stats::terms(eff))
      stopifnot(attr(mt, 'term.labels') %in% colnames(newdata))
      mf <- stats::model.frame(mt, newdata, xlev = eff$xlevels)
      stats::.checkMFClasses(attr(mt, 'dataClasses'), mf)
      mm <- stats::model.matrix(mt, mf, contrast.arg = eff$contrasts)
      stopifnot(colnames(mm) == names(coefs))
      as.numeric(mm %*% coefs)
    }
  }

  stopifnot(names(object$effect) == names(object$means$coef$effect))
  res2 <- purrr::map2(object$effect,
                      coefs_ls$effect[names(object$effect)],
                      predict_eff_i)

  ## c(res1, res2)
  if (is.null(level)) {
    purrr::reduce(c(res1, res2), `+`)
  } else {
    purrr::reduce(c(res1[level], res2), `+`)
  }
}


