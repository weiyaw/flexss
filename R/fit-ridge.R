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
#' @param K A list of numbers of interior knots for the population and subject
#'     curves. The list must contain \code{sub} and \code{pop}.
#' @param deg The degree of the spline polynomial
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
    Tmat_pop <- cbind(-1/sqrt(n_bsf_pop), stats::poly(1:n_bsf_pop, deg = deg, raw = FALSE),
                      crossprod(D_pop, solve(tcrossprod(D_pop))))
    Tmat_sub <- cbind(-1/sqrt(n_bsf_sub), stats::poly(1:n_bsf_sub, deg = deg, raw = FALSE),
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
#' @param spline A list of formulae of splines s(x, by =, K = , deg = ,
#'   ...). This routine looks up the variable from the data, before looking
#'   at the environment where the formula is defined. See `s` for more details.
#' @param random A formula ~ u1 + u2 ... specifying the random effect term.
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
  spline[[which_sub]] <- tidy_sub(spline[[which_sub]])
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
  class(fm) <- 'fsso'
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

## return a vector of indices that rearrange the data, such that data in the
##  same level of the 'by' variable are sticked together.
## fo: a formula
## data: the variable in 'fo' are lookup-ed here.
## return: a numeric vector
reorder_data <- function(fo, data) {
  qo <- match.call(s, fo[[2]])
  by <- eval(qo[['by']], envir = data, enclos = attr(fo, '.Environment'))
  order(by)
}

## detect if a formula corresponds to subject-specific curves
## fo: a formula
## return: boolean
detect_sub <- function(fo, envir = attr(fo, '.Environment')) {
  if (class(fo) == 'formula' && length(fo) == 2) {
    if (fo[[2]][[1]] == 's') {
      qo <- match.call(s, fo[[2]]) # qo = quote
    } else {
      stop("formula is not 's'.")
    }
  } else {
    stop("spline not specified as a one-sided formula.")
  }
  isTRUE(eval(qo[['is_sub']], envir = envir))
}

## return a formula for the subject-specific spline with 'block_dim' and
## 'intercept' arguments if they are not present. Also check if all the
## necessary arguments are present. By dedault, an intercept is added to subject
## curves, and block_dim is set to 'degree' (for subject curves without an
## intercept) or 'degree + 1' (for subject curves with an intercept).
## fo: a formula
## return: a formula
tidy_sub <- function(fo, envir = attr(fo, '.Environment')) { 
  qo <- fo[[2]]
  stopifnot(eval(qo[['is_sub']], envir = envir),
            !is.null(qo[['by']]),
            !is.null(qo[['knots']]),
            !is.null(qo[['degree']]))
    
  if (is.null(qo[['intercept']])) {
    qo[['intercept']] <- TRUE
    message('Add an intercept to the subject spline.')
  }
  
  if (is.null(qo[['block_dim']])) {
    ## depending on whether there's an intercept in the model
    qo[['block_dim']] <- eval(qo[['degree']], envir = envir) +
      eval(qo[['intercept']], envir = envir)
    message('block_dim of subject spline set to ', qo[['block_dim']], '.')
  }

  fo[[2]] <- qo
  fo
}

## transform formula to bs model mat (and associated details)
## fo: ~s(x, by = , knots = , deg =, ...)
## data: the data
## spl_names: name of the spline, to be append at the start of the 'by' variable
## env: envioronment which the formula should be evaluated; default to the
## environment where the formula was defined. 'x' and 'by' are evaluated
## at the data.
## return: list(model_mat, trans_mat, index, degree, type, knots)
parse_spline <- function(fo, spl_name, data, envir = attr(fo, '.Environment')) {
  if (class(fo) == 'formula' && length(fo) == 2) {
    if (fo[[2]][[1]] == 's') {
      qo <- match.call(s, fo[[2]]) # qo = quote
    } else {
      stop("formula is not 's'.")
    }
  } else {
    stop("spline must be a one-sided formula.")
  }
  
  ## ## current restriction: variables must be in the data
  ## stopifnot(as.character(qo[c('x', 'by')]) %in% c(names(data), 'NULL'))
  ## eval 'x' and 'by' at data frame first, then env
  arg_data <- purrr::map(qo[c('x', 'by')], eval, envir = data, enclos = envir)

  ## eval the rest (i.e. all except the 'x' and 'by') at env only
  rest <- setdiff(names(qo), c('', 'x', 'by'))
  arg_env <- purrr::map(qo[rest], eval, envir = envir)
  
  stopifnot(length(arg_data) == 2, names(arg_data) %in% c('', 'x', 'by'))
  c(list(call = qo), do.call(s, c(arg_data, arg_env)))
  ## s(x, by, knots, degree, intercept, penalty, type = 'bs')
}



#' Return the design matrix of a spline
#'
#' Return the design matrix of a spline.
#'
#' For non-subject-specific splines, the first 'degree + 1' terms correspond to
#' the coefficients of a polynomial of 'degree' and are not penalised. The rest
#' of the coefficients are shrunk to zero so that when all of these terms are
#' zero, the fit is a polynomial.
#'
#' Since the polynomial (which includes lines parallel to the x-axis) terms are
#' not penalised, the intercept term in the polynomial should be removed to
#' avoid spanning the same space as the (potential) main effect terms. For
#' subject-specific splines, all the coefficients are penalised, so including an
#' intercept will not result in an unidentifiable model.
#'
#' Due to the way which B-spline basis is parameterised and transformed, it is
#' advisable to treat the first 'degree + 1' coefficients of the
#' subject-specific spline as correlated Gaussian random effects (i.e. with full
#' covariance matrix). The rest of the spline coefficients (subject-specific or
#' otherwise) are independent Gaussian. 'block_dim' is only applicable for the
#' subject-specific spline.
#'
#' @param x a variable of the design points.
#' @param by a grouping variable, or if NULL, a common spline is fitted.
#' @param knots the number of internal knots.
#' @param degree the degree of spline.
#' @param intercept TRUE/FALSE whether an intercept should be included in the
#'   spline.
#' @param type the type of basis function. Currently only B-spline ('bs') is
#'   supported.
#' @param is_sub is this the subject-specific curves?
#' @param block_dim the (first) number of terms in the subject-specific curves
#'   that should be treated as a block, rather than i.i.d.
s <- function(x, by = NULL, knots = NULL, degree = NULL, intercept = NULL,
              type = 'bs', is_sub = FALSE, block_dim = NULL) {

  res <- list(degree = degree,
              intercept = intercept,
              type = type,
              is_sub = is_sub,
              block_dim = block_dim)

  ## a list of positions of the data points for each group
  if (type != 'bs') {
    stop('basis type not implemented yet.')
  }

  if (!(is.null(by) || is.factor(by))) {
    stop('by variable must be a factor.')
  }

  stopifnot(is.numeric(knots), length(knots) == 1)
  dmatinfo <- get_design_bs(x, knots, degree)

  ## if intercept is true, add an intercept
  if (!is.null(intercept) && intercept) {
    res$trans_mat <- get_transform_bs(knots + degree + 1, degree + 1)
    cp <- cbind(matrix(0, knots, degree + 1), diag(knots)) # compact penalty
  } else {
    ## remove the intercept of the spline
    res$trans_mat <- get_transform_bs(knots + degree + 1, degree + 1)[, -1]
    cp <- cbind(matrix(0, knots, degree), diag(knots)) # compact penalty
  }

  cmm <- dmatinfo$design %*% res$trans_mat         # compact model matrix
  stopifnot(is.null(by) || length(by) == NROW(cmm)) # otherwise will screw up predict(.)
  if (is_sub) {
    ## check if sub index is continuous, bcoz assuming Bmat$sub is block diagonal
    index <- split(1:length(by), by)
    level <- names(index)
    stopifnot(purrr::map_lgl(index,
                             ~all(.x == seq.int(.x[1], len = length(.x)))))
    res$model_mat <- purrr::map(index, ~cmm[.x, ])
    attr(res$model_mat, 'block_dim') <- block_dim
    attr(res$model_mat, 'index') <- index
  } else {
    if (is.null(by)) {
      level <- 'mean'
      res$model_mat <- cmm
    } else {
      level <- unique(by)
      smm_ls <- purrr::map(level, # sparse model matrix
                           ~Matrix::Matrix(cmm * (by == .x), sparse = TRUE))
      res$model_mat <- do.call(cbind, smm_ls)
    }
    attr(res$model_mat, 'penalty') <- cp
    stopifnot(NCOL(cmm) == NCOL(cp))
  }
  res$knots <- dmatinfo$knots

  ## these four attributes are required for Gibbs sampler
  attr(res$model_mat, 'spl_dim') <- NCOL(cmm)
  attr(res$model_mat, 'is_sub') <- is_sub
  attr(res$model_mat, 'level') <- as.character(level)
  res
}

## parse the formula to return a list of a term object and a vector of response.
parse_response <- function(fo, data, envir = attr(fo, '.Environment')) {
  if (is.null(fo)) {
    stop("response required.")
  }
  mt <- stats::terms(fo, data = data)[0] # mt: model terms
  if (attr(mt, 'response') == 0) {
    stop("response required.")
  }
  
  ## ## current restriction: variables must be in the data
  ## stopifnot(as.character(attr(mt, 'variables')[[2]]) %in% names(data))

  list(terms = mt,
       y = eval(attr(mt, 'variables')[[2]], data, envir))
}

## parse the formula to return a list of a term object, model matrix, and some
## attributes useful for predicting new data points.
parse_effect <- function(fo, data, envir = attr(fo, '.Environment')) {
  if (is.null(fo)) {
    NULL
  } else {
    ## mt: model term
    mt <- stats::delete.response(stats::terms(fo, data = data))
    stopifnot(attr(mt, 'response') == 0) # must not have response term in the formula
    
    eval_variable <- eval(attr(mt, 'variable'), data[1, ], envir)
    attr(mt, 'dataClasses') <- purrr::map_chr(eval_variable, stats::.MFclass)
    names(attr(mt, 'dataClasses')) <- as.character(attr(mt, 'variable')[-1])

    ## ## current restriction: variables must be in the data
    ## stopifnot(attr(mt, 'term.labels') %in% names(data))

    ## mm: model matrix, all terms must exist in data
    mm <- stats::model.matrix(mt, data)

    list(terms = mt,
         model_mat = mm,
         assign = attr(mm, 'assign'),
         contrasts = attr(mm, 'contrasts'),
         xlevels = stats::.getXlevels(mt, data))
  }
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
#' @param object an 'fsso' object.
#' @param newdata a data frame with all the covariate required for prediction.
#' @param level a string specify the name of the spline. Only the component
#'   correponding to that spline is returned together with the effect terms).
#' @param fun a function applying to the samples of each marginal posterior of
#'   the regression coefs. By default, it's `mean`.
#' 
#' @export
predict.fsso <- function(object, newdata = NULL, level = NULL, fun = NULL) {
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


