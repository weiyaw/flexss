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



#' Return the design matrix (and other objects) for fitting a penalised spline
#'
#' Return the design matrix of a spline, and objects such as the penalty matrix
#' that are required for running the spline fitting engine
#' (bayes_ridge_semi). For the time being, the only supported basis is
#' 'B-spline'. The design matrix is transformed such that the model can be
#' fitted as a linear mixed-effect models.
#'
#' The 'by' argument specifies the factor variable that interacts with the 'x'
#' covariate. Basically, this means that whether we should have a spline for
#' each level of the factor. If there should not be any interaction, we can
#' specify a dummy factor variable, e.g. as.factor(1).
#' 
#' This function is unlikely to be called by the users directly but rather by
#' the fitting routine.
#'
#' @details # Working assumption
#'
#' ## Non-subject-specific splines
#'   The polynomial coeffcients have a flat prior, whereas the rest have a
#'   zero-mean independent Gaussian prior. This implies that the splines are
#'   shrunk to the polynomial of degree `degree`. Since the polynomial space
#'   (which includes an intercept) spanned by a spline are not penalised, the
#'   intercept term is removed to avoid spanning the same space as the
#'   (possible) main effect terms and resulting in an unidentifiable model. This
#'   usually means that the users should fit the intercepts as fixed main
#'   effects. For example, if you have `s(x1, by = b1)`, then you should have
#'   something like `y ~ b1` as you fixed (main effect) term.
#'
#' ## Subject-specific splines
#'   All regression coefficients are penalised
#'   towards 0, i.e. has a zero-mean Gaussian prior. Intercept is included. Due
#'   to the LMM parameterisation of B-spline basis, the first 'degree + 1'
#'   coefficients of the subject-specific spline are Gaussian with a full
#'   covariance matrix. The rest of the coefficients are independent Gaussian.
#' 
#'
#' @param x a numeric vector of the design points.
#' @param by a grouping variable and must be a factor. Vectors of length one are
#'   expanded to `length(x)`. Use a dummy variable, e.g. as.factor(1), if there
#'   is only one group.
#' @param knots the number of internal knots.
#' @param degree the degree of spline.
#' @param type the type of basis function. Currently only B-spline ('bs') is
#'   supported.
#' @param is_sub is this the subject-specific curves?

s <- function(x, by = as.factor(1), knots = min(max(1, length(x)/4), 35), degree = 2,
              type = 'bs', is_sub = FALSE) {

  res <- list(degree = degree,
              type = type,
              is_sub = is_sub)

  ## a list of positions of the data points for each group
  if (type != 'bs') {
    stop('basis type not implemented yet.')
  }

  if (!is.factor(by)) {
    stop('by variable must be a factor.')
  }

  ## If there is only 1 population
  if (!is_sub && length(by) == 1) {
    by <- rep(by, times = length(x))
  } else if (length(x) != length(by)) {
    stop("length of 'x' and 'by' should match.")
  }
  stopifnot(is.numeric(knots), length(knots) == 1)
  dmatinfo <- get_design_bs(x, knots, degree)
  res$knots <- dmatinfo$knots

  if (is_sub) {
    res$trans_mat <- get_transform_bs(knots + degree + 1, degree + 1)
    cmm <- dmatinfo$design %*% res$trans_mat # compact model matrix

    ## check if sub index is continuous, bcoz assuming Bmat$sub is block diagonal
    index <- split(1:length(by), by)
    level <- names(index)
    stopifnot(purrr::map_lgl(index,
                             ~all(.x == seq.int(.x[1], len = length(.x)))))
    res$model_mat <- purrr::map(index, ~cmm[.x, ])

    attr(res$model_mat, 'block_dim') <- degree + 1
    message('block_dim of subject spline set to ', attr(res$model_mat, 'block_dim'), '.')
    attr(res$model_mat, 'index') <- index
  } else {

    ## Remove the intercept of the spline to ensure model identifiability. The
    ## intercept shall be fitted as a fixed main effect.
    res$trans_mat <- get_transform_bs(knots + degree + 1, degree + 1)[, -1]
    cmm <- dmatinfo$design %*% res$trans_mat           # compact model matrix
    cp <- cbind(matrix(0, knots, degree), diag(knots)) # compact penalty
    level <- levels(by) # this has to match columns in model_mat
    if (length(level) == 1) {
      res$model_mat <- cmm
    } else {
      smm_ls <- purrr::map(level, # sparse model matrix
                           ~Matrix::Matrix(cmm * (by == .x), sparse = TRUE))
      res$model_mat <- do.call(cbind, smm_ls)
    }
    attr(res$model_mat, 'penalty') <- cp
    stopifnot(NCOL(cmm) == NCOL(cp)) # sanity check
  }

  ## sanity check
  stopifnot(length(by) == NROW(cmm)) # otherwise will screw up predict(.)

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

