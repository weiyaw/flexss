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
## necessary arguments are present. By default, an intercept is added to subject
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

  ## these three attributes are required for Gibbs sampler
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

