## inverse of a symmetrical, positive definite matrix
pdinv <- function(x, tol = .Machine$double.eps * 100000) {
  ## HERE: POTENTIAL NUMERICAL ISSUE, MIGHT NEED TO ADJUST TOLERANCE
  if (!Matrix::isSymmetric(x, tol = tol)) {
    if (any((x - Matrix::t(x)) > tol)) {
      message('asymmetrical pd matrix.')
      x <- (x + Matrix::t(x)) / 2      
    }
  }
  ## use 'dpoMatrix' only if we can guarantee x is pd
  ## xpd <- as(Matrix::forceSymmetric(x), 'dpoMatrix')
  Matrix::forceSymmetric(Matrix::solve(x))
}

## product of the quadratic, x %*% A %*% t(x)
tquadprod <- function(x, A) {
  res <- Matrix::tcrossprod(x %*% A, x)
  ## stopifnot(Matrix::isSymmetric(res, tol = .Machine$double.eps * 10000))
  Matrix::forceSymmetric(res)
}

## is this a precision matrix (in the context of this model)? The precision can
## be semi-positive definite.
is_precision <- function(x) {
  all(isSymmetric(x), eigen(x, TRUE)$values >= 0)
}

check_prec_v4 <- function(prec, para) {
  ## display specified precision
  if (!is.null(prec$theta)) {
    message('prec$theta specified.')
    stopifnot(is.matrix(prec$theta),
              dim(prec$theta) == para$dim_theta,
              is_precision(prec$theta))
  }

  if (!is.null(prec$delta)) {
    message('prec$delta specified.')
    stopifnot(is.matrix(prec$delta),
              dim(prec$delta) == para$dim_delta,
              is_precision(prec$delta))
  }

  if (!is.null(prec$beta)) {
    if (is.null(para$dim_beta)) {
      message("prec$beta specified but irrelavent. No fixed/random effects.")
    } else {
      message('prec$beta specified.')
      stopifnot(is.matrix(prec$beta),
                dim(prec$beta) == para$dim_beta,
                is_precision(prec$beta))
    }
  }

  if (!is.null(prec$eps)) {
    message('prec$eps specified.')
    stopifnot(is.numeric(prec$eps),
              length(prec$eps) == 1,
              prec$eps > 0)
  }
}


## initialise coef with some random Gaussian rv.
initialise_v4 <- function(Bmat, Xmat, y) {

  ## penalised least squares with lambda = 100, and add some gaussian noise
  pls <- function(x) {
    if (is.list(x) && attr(x, 'is_sub')) {
      beta_pls <- 0
    } else if (is.matrix(x)) {
      beta_pls <- solve(crossprod(x) + 100 * crossprod(attr(x, 'penalty'))) %*%
        crossprod(x, y)
    } else {
      stop('Cannot initialise.')
    }
    
    ## add some noise into the pls
    beta_pls + 
      matrix(stats::rnorm(attr(x, 'spl_dim') * length(attr(x, 'level'))),
             nrow = attr(x, 'spl_dim'),
             ncol = length(attr(x, 'level')))
    
  }

  res <- list(spl = purrr::map(Bmat, ~structure(pls(.x),
                                                is_sub = attr(.x, 'is_sub'),
                                                penalty = attr(.x, 'penalty'),
                                                block_dim = attr(Bmat[[1]], 'block_dim'))),
              eff = purrr::map(Xmat, ~structure(stats::rnorm(NCOL(.x)))))

  res$yhat <- y + stats::rnorm(y)
  res
}

initialise_prec_v4 <- function(Bmat, Xmat, y) {
  spl_prec <- function(x) {
    if (attr(x, 'is_sub') & is.null(attr(x, 'penalty'))) {
      diag(attr(x, 'spl_dim'))
    } else if (!attr(x, 'is_sub') & !is.null(attr(x, 'penalty'))) {
      stopifnot(NCOL(attr(x, 'penalty')) == attr(x, 'spl_dim'))
      crossprod(attr(x, 'penalty')) * 100
    } else {
      stop('Cannot initialise.')
    }
  }
  
  eff_prec <- function(x) {
    diag(NCOL(x))
  }

  list(spl = purrr::map(Bmat, spl_prec),
       eff = purrr::map(Xmat, eff_prec),
       eps = 1 / stats::sd(y))
}


## Bmat: list of matrix
## return: list of matrix
calc_xBs_v4 <- function(Bmat_sub) {
  purrr::map(Bmat_sub, Matrix::crossprod)
}

## Phi: matrix / list of matrix if Phi_1
## Xmat: matrix
## return: matrix
calc_PhixX_v4 <- function(Phi, Xmat) {
  if (is.list(Phi)) {
    Phi <- Matrix::.bdiag(Phi)
  }
  Phi %*% Xmat
}

## Xmat: matrix
## PhixX: matrix
## kprec: list(beta: matrix, eps: numeric)
## return: matrix
calc_Qmat_inv_v4 <- function(Xmat, PhixX, Sigma_inv, eps_inv) {
  if (is.null(Xmat) || is.null(PhixX)) {
    NULL
  } else {
    XtX <- Matrix::crossprod(Xmat, PhixX)
    sigma_ratio <- sign(Sigma_inv) * exp(log(abs(Sigma_inv)) - log(eps_inv))
    
    spl_dim <- attr(Xmat, 'spl_dim')
    if (is.null(spl_dim)) {
      ## if it's an effect Xmat (i.e. Xmat in the paper)
      ## if it's a fixed effect, sigma_ratio should be 0
      ## eps_inv is strictly positive
      pdinv(XtX + sigma_ratio)
    } else {
      ## if it's a spline Xmat (i.e. Bmat in the paper)
      for (i in seq_along(attr(Xmat, 'level'))) {
        idx <- seq.int(to = spl_dim * i, length.out = spl_dim)
        XtX[idx, idx] <- XtX[idx, idx] + sigma_ratio
      }
      pdinv(XtX)
      ## }
    }
  }
}

## assume: when computing Phi2, Phi1 is a list but converted to a bdiag matrix
## (i.e. Bmat_sub assumed to be block diagonal)
## Phi: matrix / list of matrix if Phi_2
## PhixB: matrix / list of matrix if Phi_1
## Qmat_inv: matrix / list of matrix if Phi_1
## return: matrix / list of matrix if Phi_1
calc_Phi_v4 <- function(Phi, PhixX, Qmat_inv) {
  if (is.null(Phi)) {
    ## calc base case Phi_1, PhixX and Qmat_inv are lists
    stopifnot(is.list(PhixX), is.list(Qmat_inv))
    calc_Phi0i <- function(PhixX, Qmat_inv) {
      res <- -1 * tquadprod(PhixX, Qmat_inv)
      res[row(res) == col(res)] <- Matrix::diag(res, names = FALSE) + 1
      ## stopifnot(Matrix::isSymmetric(res, tol = .Machine$double.eps * 10000))
      Matrix::forceSymmetric(res)
    }
    purrr::map2(PhixX, Qmat_inv, calc_Phi0i)
  } else {
    if (is.list(Phi)) {
      ## when computing Phi_2, Phi_1 is a list
      Phi <- Matrix::.bdiag(Phi)
    }
    ## computing Phi_2, 3, 4, ...
    Phi - tquadprod(PhixX, Qmat_inv)
  }
}

## calculate partial residual
## presid: col matrix
## Xb: col matrix
calc_presid_v4 <- function(presid, Xb) {
  if (is.null(Xb)) {
    presid
  } else {
    ## if (is.list(presid) && is.list(Xb)) {
    ##   stopifnot(is.list(presid), is.list(Xb))
    ##   purrr::map2(presid, Xb, ~as.matrix(.x - .y))
    ## } else {
    presid - Xb
    ## }
  }
}

## calculate Rmat
## presid: col matrix
## PhixX: matrix / list of matrix (i.e. bdiag mat)
## return: row matrix / list of row matrix (if bdiag mat)
calc_Rmat_v4 <- function(presid, PhixX) {
  if (is.list(PhixX)) {
    ## in this scenario, PhixX is just Bmat of subject terms
    n <- vapply(PhixX, NROW, 1L)
    stopifnot(sum(n) == length(presid), length(n) == length(PhixX))
    by <- rep(seq_along(PhixX), times = n)
    ## split presid into a list of numeric vectors
    presid_ls <- split(as.numeric(presid), by)
    ## calculate Rmat for each subject
    purrr::map2(presid_ls, PhixX, Matrix::crossprod)
  } else {
    Matrix::crossprod(presid, PhixX)
  }
}

## draw samples from joint conditional posterior
## Qmat_inv: matrix / list of matrix
## Rmat: row matrix / list of row matrix
## eps_inv: numeric
## return: numeric / list of numeric
gen_coefs <- function(Qmat_inv, Rmat, eps_inv) {
  gen_coefs_i <- function(x, y) {
    mu <- as.numeric(Matrix::tcrossprod(x, y))
    ## eps_inv is strictly positive
    Sigma <- as.matrix(sign(x) * exp(log(abs(x)) - log(eps_inv)))
    MASS::mvrnorm(1, mu, Sigma)
  }
  if (is.list(Qmat_inv) || is.list(Rmat)) {
    stopifnot(is.list(Qmat_inv), is.list(Rmat))
    purrr::map2(Qmat_inv, Rmat, gen_coefs_i)
  } else {
    gen_coefs_i(Qmat_inv, Rmat)
  }
}

## Xmat: matrix / list of matrix
## beta: col matrix / list of col matrix
## return: col matrix / list of matrix
calc_Xb_v4 <- function(Xmat, beta) {
  stopifnot(!is.null(Xmat), !is.null(beta))
  if (is.list(Xmat) || is.list(beta)) {
    stopifnot(is.list(Xmat), is.list(beta))
    purrr::map2(Xmat, beta, ~as.numeric(.x %*% .y))
  } else {
    Xmat %*% beta
  }
}

update_yhat_v4 <- function(f, g, Xb) {
  stopifnot(length(f) == length(g))
  if (is.null(Xb)) {
    f + g
  } else {
    stopifnot(length(f) == length(Xb))
    f + g + Xb
  }
}


## datals: list(Bmat_pop: matrix, Bmat_sub: list of matrix, Xmat: matrix, y: col matrix)
## y: col matrix
## Bmat: list of matrix. term for subject curves must be at the front.
## Xmat: list of matrix
## xBs: list of matrix (crossprod(Bmat$sub) for each i)
## kprec: list(spl: list of matrix, eff: list of matrix, eps: numeric)

#' Update coefficients
#'
#' @param y An n by 1 column matrix of the responce
#' @param Bmat A list of spline term design matrices. The term for subject
#'   curves must be at the front. Required attributes of each matrix: spl_dim
#'   (num), level (char), is_sub (bool), penalty (matrix, non-sub spline only),
#'   block_dim (num, sub spline only)
#' @param Xmat
#' 
#' @return precision matrix for the coefs
update_coefs_v4 <- function(y, Bmat, Xmat, xBs, kprec, debug = FALSE) {

  ## term for subject curves must be at the front.
  stopifnot(attr(Bmat[[1]], 'is_sub'))
  
  allmat <- c(Bmat, Xmat)
  allprec <- c(kprec$spl, kprec$eff)
  stopifnot(names(allmat) == names(allprec))

  ## add 's' and 'e' at the front to prevent name crash
  allnames <- c(paste0('s', names(Bmat)), paste0('e', names(Xmat)))
  stopifnot(length(allnames) == length(unique(allnames)))
  names(allmat) <- allnames
  names(allprec) <- allnames

  Qmat_inv <- purrr::map(allmat, ~NULL)
  PhixX <- purrr::map(allmat, ~NULL)
  Rmat <- purrr::map(allmat, ~NULL)
  coefs <- purrr::map(allmat, ~NULL)
  
  ## subject-specific terms
  PhixX[[1]] <- Bmat[[1]]
  sigma_sub_ratio <- sign(allprec[[1]]) * exp(log(abs(allprec[[1]])) - log(kprec$eps))
  Qmat_inv[[1]] <- purrr::map(xBs, ~pdinv(.x + sigma_sub_ratio))
  Phi <- NULL

  ## calculate Qmat_inv
  for (l in seq_along(allnames)[-1]) {     
    ct <- allnames[l] # ct: current term
    pt <- allnames[l-1] # pt: previous term
    Phi <- calc_Phi_v4(Phi, PhixX[[pt]], Qmat_inv[[pt]])
    PhixX[[ct]] <- calc_PhixX_v4(Phi, allmat[[ct]])
    Qmat_inv[[ct]] <- calc_Qmat_inv_v4(allmat[[ct]], PhixX[[ct]], allprec[[ct]], kprec$eps)
  }
  
  presid <- y
  Xb <- NULL
  ## calculate Rmat
  for (ct in rev(allnames)) {
    presid <- calc_presid_v4(presid, Xb)
    Rmat[[ct]] <- calc_Rmat_v4(presid, PhixX[[ct]])    
    coefs[[ct]] <- gen_coefs(Qmat_inv[[ct]], Rmat[[ct]], kprec$eps)
    Xb <- calc_Xb_v4(allmat[[ct]], coefs[[ct]])
  }
  
  res <- list(spl = coefs[grepl('^s', names(coefs))],
              eff = coefs[grepl('^e', names(coefs))])

  res <- purrr::map(res, ~setNames(.x, sub('^(s|e)', '', names(.x))))
  ## here, Xb is the linear prediction from the subject-specifc splines
  res$yhat <- y - presid + unlist(Xb, use.names = FALSE)
  
  ## convert spl coefs to matrices
  for (l in names(res$spl)) {
    if (attr(Bmat[[l]], 'is_sub')) {
      # subject splines
      res$spl[[l]] <- structure(do.call(cbind, res$spl[[l]]),
                                is_sub = TRUE,
                                block_dim = attr(Bmat[[l]], 'block_dim'))
    } else {
      ## return penalty for generic spline
      res$spl[[l]] <- structure(matrix(res$spl[[l]],
                                       nrow = attr(Bmat[[l]], 'spl_dim'),
                                       dimnames = list(NULL, attr(Bmat[[l]], 'level'))),
                                is_sub = FALSE,
                                penalty = attr(Bmat[[l]], 'penalty'))
    }
  }
  res
}

get_coef_container_v4 <- function(para, size) {
  
  spl_array <- purrr::map2(para$spl_dims, para$spl_level, 
                           ~array(NA, c(.x, length(.y), size),
                                  dimnames = list(NULL, .y, NULL)))

  eff_matrix <- purrr::map2(para$eff_dims, para$eff_level,
                            ~matrix(NA, .x, size, dimnames = list(.y, NULL)))
  
  list(spline = spl_array, effect = eff_matrix)
}

get_prec_container_v4 <- function(para, size) {
  spl_array <- purrr::map(para$spl_dims, ~array(NA, c(.x, .x, size)))
  eff_array <- purrr::map(para$eff_dims, ~array(NA, c(.x, .x, size)))
  eps_vector <- rep(NA, size)

  list(spline = spl_array, effect = eff_array, eps = eps_vector)
}


#' Get dims and names of the spline and effect terms.
#'
#' @param Bmat,Xmat Design matrices of the spline and effect terms.
#' 
#' @return A list of dims and levels
get_para_v4 <- function(Bmat, Xmat) {
  list(spl_dims = purrr::map(Bmat, ~attr(.x, 'spl_dim')),
       eff_dims = purrr::map(Xmat, NCOL),
       spl_level = purrr::map(Bmat, ~attr(.x, 'level')),
       eff_level = purrr::map(Xmat, colnames))
}

## x: a vector or matrix, to be squared and sum
update_with_gamma <- function(x, a, b) {
  stopifnot(!is.null(x), !is.null(a), !is.null(b))
  shape <- 0.5 * length(x) + a
  rate <- 0.5 * sum(x^2) + b
  stats::rgamma(1, shape = shape, rate = rate)
}

## x: a matrix of col vectors x, to be tcrossprod
update_with_wishart <- function(x, v, lambda) {
  stopifnot(!is.null(x), !is.null(v), !is.null(lambda),
            NCOL(lambda) == NROW(x))
  df <- v + NCOL(x)
  scale <- lambda + Matrix::tcrossprod(x)
  inv_scale <- as.matrix(pdinv(scale))
  stats::rWishart(1, df = df, Sigma = inv_scale)[, , 1]
}

#' Update precision of the spline terms
#'
#' @param coefs Matrix, each column the coef of (a subject|a level in a
#'   factor). `attr(coefs, 'penalty')`: a full-row rank penalty matrix (for
#'   ordinary splines). Kmat %*% coefs are penalised (i.e. are
#'   random). `attr(coefs, 'block_dim')`: the dimension corresponding to the
#'   block cov matrix (for subject splines)
#' 
#' @return precision matrix for the coefs
update_prec_spline <- function(coefs, prior) {
  if (is.null(prior)) {
    NULL
  } else {
    if (attr(coefs, 'is_sub')) {
      block_dim <- attr(coefs, 'block_dim')
      ## when coefs are for subject splines
      if (block_dim > NROW(coefs) || block_dim < 0) {
        stop('block_dim out of bound.')
      }
      block <- NA
      iid <- NA
      ## block precision term
      if (block_dim > 0) {
        coefs_block <- coefs[seq(1, block_dim), , drop = FALSE]
        block <- update_with_wishart(coefs_block, prior$v, prior$lambda)
      }
      ## iid precision term
      if (block_dim < NROW(coefs)) {
        coefs_iid <- coefs[seq(block_dim + 1, NROW(coefs)), , drop = FALSE]
        iid <- update_with_gamma(coefs_iid, prior$a, prior$b)
      }
      
      stopifnot(length(iid) == 1)        # assuming iid is a scaler
      `[<-`(diag(iid, NROW(coefs)), seq_len(block_dim), seq_len(block_dim), block)
    } else {
      ## when coefs are for ordinary splines
      Kmat <- attr(coefs, 'penalty')
      stopifnot(!is.null(Kmat))
      x <- Kmat %*% coefs
      update_with_gamma(x, prior$a, prior$b) * Matrix::crossprod(Kmat)
    }
  }
}
## yhat, y: list of col vectors/vectors with the same names
## return precision of the residual
update_prec_eps <- function(yhat, y, prior) {
  stopifnot(names(y) == names(yhat))
  resids <- y - yhat
  ## resids <- unlist(y) - unlist(yhat)
  update_with_gamma(resids, prior$a, prior$b)
}

## coefs: col vector/vector
## prior: list of prior hyperparameter
## if prior is NULL, return a zero matrix. (i.e. coefs are fixed effects)
update_prec_effect <- function(coefs, prior) {
  stopifnot(NCOL(coefs) == 1)
  if (is.null(prior)) {
    matrix(0, NROW(coefs), NROW(coefs))
  } else {
    prec <- update_with_gamma(coefs, prior$a, prior$b)
    diag(prec, NROW(coefs))
  }
}

update_precs_v4 <- function(kcoef, y, prior_ls, init_prec = NULL) {
  if (is.null(init_prec)) {
    stopifnot(names(kcoef$spl) == names(prior_ls$spl),
              names(kcoef$eff) == names(prior_ls$eff))    
    list(spl = purrr::map2(kcoef$spl, prior_ls$spl, update_prec_spline),
         eff = purrr::map2(kcoef$eff, prior_ls$eff, update_prec_effect),
         eps = update_prec_eps(kcoef$yhat, y, prior_ls$eps))
  } else {
    init_prec
  }
}

recurse <- function(x, f) {
  if (is.list(x)) {
    lapply(x, recurse, f = f)
  } else {
    f(x)
  }
}
  
#' Check list structure
#'
#' Check if all the lists have the same structure.
#'
#' @param ... Lists to be checked.
#'
#' @return Boolean
is_match_list <- function(...) {
  ## prune the leaves of lists, leaving only structures
  struct <- lapply(list(...), recurse, f = function(x) NA)
  ## check if the structures are the same
  all(vapply(struct, identical, TRUE, y = struct[[1]]))
}

#' Get posterior statistics
#'
#' Calculate posterior statistics from an array or vector of samples. The samples are
#' populated along the last dimension, e.g. columns of a matrix are the samples;
#' the depth of a 3D array are the samples.
#'
#' @param samples A vector, matrix or array of samples
#' @param fun A function for statitics calculation. Support `tidy` style of
#'   function specification.
#'
#' @return A numeric value, vector or matrix
#' 
#' @name pstats
NULL

#' @rdname pstats
#' 
#' @details `pstats` calculates posterior statistics given by `fun`.
pstats_v4 <- function(samples, fun) {
  fun <- purrr::as_mapper(fun)
  if (is.array(samples)) {
    apply(samples, seq(1, length(dim(samples)) - 1), fun)
  } else if (is.vector(samples, mode = 'numeric')) {
    fun(samples)
  ## } else if (is.list(samples)) {
  ##   purrr::map(samples, pstats_v4, fun = fun)
  } else {
    stop("Invalid samples structure.")
  }
}

#' @rdname pstats
#' 
#' @details `pmean` calculates posterior means, and faster than pstats(samples, mean).
pmean_v4 <- function(samples) {
  if (is.array(samples)) {
    rowMeans(samples, dims = length(dim(samples)) - 1)
  } else if (is.vector(samples, mode = 'numeric')) {
    mean(samples)
  ## } else if (is.list(samples)) {
  ##   purrr::map(samples, pmean_v4)
  } else {
    stop("Invalid samples structure.")
  }
}



#' Convert precision to variance/covariance
#'
#' Convert precision (matrices) to variance (covariance matrices).
#'
#' @param prec a vector of precision, an array of precision matrices
#'
#' @return variance or covariance matrix
prec_to_cov <- function(prec) {
  if (is.vector(prec, mode = 'numeric')) {
    1 / prec
  } else if (is.array(prec) && length(dim(prec)) == 3) {
    size <- dim(prec)[3]
    for (i in 1:size) {
      prec[, , i] <- chol2inv(chol(prec[, , i]))
    }
    prec
  ## } else if (is.list(prec)) {
  ##   purrr::map(prec, prec_to_cov)
  } else {
    stop("Invalid precision structure.")
  }
}


check_Bmat <- function(Bmat) {
  
  check_each_Bmat <- function(x, x_name) {
    if (is.null(attr(x, 'spl_dim'))) {
      stop('spl_dim not specified in ', x_name, '.')
    }
    if (is.null(attr(x, 'is_sub'))) {
      stop('is_sub not specified in ', x_name, '.')
    } 
    if (is.null(attr(x, 'level'))) {
      stop('level not specified in ', x_name, '.')
    } 
    if (!is.character(attr(x, 'level'))) {
      ## This is to ensure that numeric levels will not mess up vector indexing.
      stop('level in ', x_name, ' has to be a character (vector).')
    } 

    if (attr(x, 'is_sub')) {
      if (!(is.list(x) &&
              all(purrr::map_lgl(x, ~is.matrix(.x) || methods::is(.x, 'Matrix'))))) {
        stop(x_name, ' must be a list of matrices.')
      }
      if (!all(purrr::map_dbl(x, NCOL) == attr(x, 'spl_dim'))) {
        stop('number of cols of ', x_name, ' mismatch with spl_dim and level.')
      } 
      if (!identical(sort(names(x)), sort(attr(x, 'level')))) {
        stop('names of the level in the subject term and list of matrices must match up.')
      }

      if (is.null(attr(x, 'block_dim'))) {
        stop('missing block_dim in the ', x_name, '.')
      }
    } else {
      if (!(is.matrix(x) || methods::is(x, 'Matrix'))) {
        stop(x_name, ' must be a matrix.')
      }
      if (NCOL(x) != (attr(x, 'spl_dim') * length(attr(x, 'level')))) {
        stop('incompatible column dimension for ', x_name, '.')
      } 
      if (is.null(attr(x, 'penalty'))) {
        stop('penalty not specified in ', x_name, '.')
      } else {
        if (abs(Matrix::det(Matrix::tcrossprod(attr(x, 'penalty')))) < 1e-10) {
          warning('rows of Kmat may not be independent.')
        }
      }
    }
  }

  if (any(purrr::map_lgl(names(Bmat), is.null))) {
    stop('all entries of Bmat must have a name.')
  }
  if (sum(purrr::map_lgl(Bmat, ~attr(.x, 'is_sub'))) != 1) {
    stop('must include one and only one term for subject-specific curves.')
  }
  purrr::iwalk(Bmat, check_each_Bmat)
}

check_Xmat <- function(Xmat) {

  if (any(purrr::map_lgl(names(Xmat), is.null))) {
    stop('all entries of Xmat must have a name.')
  }
}


## This is a Gibbs sampler v3 for longitudinal Bayesian semiparametric ridge.
## Update two block of parameters: variance and coefs

## Requirements: y, grp, Bmat (df), Xmat (df), Kmat, block_dim, ranef

## Algorithms paremeters: burn, size
## Extras: init (list of matrices with 'pop' and 'sub'), prior (see 'check_prior')


#' Bayesian ridge for longitudinal semiparemetric models
#'
#' This is a Gibbs sampler v4 for fitting Bayesian longitudinal semiparametric
#' models. It assumes that there are multiple populations in the model and each
#' population is modelled by a population 'mean' curve. Within each population,
#' they are multiple subjects, and each subject are modelled by a 'subject'
#' curve. The subject curves are treated as deviations from their respective
#' mean curves. On top of that, fixed or random effects can be added to the
#' model. This is useful when, for example, a particular treatment was applied
#' to some of the populations or subjects, and the user is interested in the
#' effect of the treatment.
#'
#' THIS FUNCTION IS FOR INTERNAL USE ONLY AND NEVER MEANT TO BE EXPORTED.
#' 
#' The mean and subject curves are modelled as linear (in statistical sence, not
#' only stright lines). This includes polynomials and splines. The attributes of
#' their design matrices (i.e. Bmat) are
#'
#' 'spl_dim' is the dimension of the splines. This is actually redundant (can be
#' inferred from 'level' and the dims of Bmat), but is still included for
#' verification and for easy query of the spline dimension.
#' 
#' 'is_sub' is a boolean that specifies whether the spline term is a
#' subject-specific deviations. One of the spline must be a subject spline.
#'
#' 'level' is the name of each population or subject. The order of the names
#' should match the corresponding column entries in Bmat/list of Bmat. Use a
#' dummy name if there is only one level.
#'
#' 'block_dim' is the dimension of the subject deviations where the precision
#' should be considered as a block. See paper for more details.
#' 
#' 'penalty' relates to the precision of the prior of the mean curve
#'   coefficient, which is a multiple of crossprod('penalty'), i.e. crossprod('penalty'
#'   %*% \theta_i) is shrunk to 0.
#'
#' The semiparametric model is
#'
#' y = Bmat[[1]] %*% \theta_1 + ... + Xmat[[1]] %*% \beta_1 + ... + Bmat[[sub]]
#' %*% \delta + \epsilon
#'
#' y is a vector of observation, Bmat[[...]] %*% \theta are the splines that
#' constitute the mean curve, Xmat[[...]] %*% \beta are the fixed/random
#' effects, Bmat[[sub]] %*% \deltais the subject-specific deviations, and
#' \epsilon are Gaussian errors.
#'
#' The coefficients \theta, \delta, \beta and \epsilon are all Gaussian, but
#' their covariance structures are not necessarily diagonal. Consult the
#' manual/paper for more info.
#'
#' The sampler first updates the coefficients, then the precisions. The joint
#' conditional posterior of all the coefficients is often highly correlated, but
#' fortunately a closed-form expression is available and is utilised here. This
#' implies that the convergence of this Gibbs sampler is quick and does not
#' require burning many samples in the beginning. The starting values of the
#' sampler are precisions with a reasonably large value (hence a large penalty).
#' 
#' @param y A vector of the response vector.
#' @param Bmat A named list of design matrices (or a list of design matrices for
#'   each subject for the subject spline) for the splines terms with at least
#'   the following attributes: 'spl_dim', 'is_sub', 'level'. For
#'   subject-specific splines (i.e. is_sub = TRUE), it is a block-diagonal
#'   matrix and specified as a list of matrices corresponding to each diagonal
#'   element. The list should also have a 'block_dim' attribute. For other
#'   splines, it is a matrix and should have a 'penalty' attribute. At least an
#'   element of Bmat must be 'is_sub = TRUE'. See details.
#' @param Xmat A named list of design matrices for the non-spline terms,
#'   e.g. fixed and random effects.
#' @param prior unknown yet
#' @param prec An optional named list with three elements 'spl', 'eff' and
#'   'eps'. If supplied, the precisions of the model are fixed. Each of the
#'   'spl' and 'eff' terms is a named list of precisions that correspond to
#'   'Bmat' and 'Xmat' respectively, and the names of the list must match. All
#'   the elements are square matrices with the dimension matching 'spl_dim' of
#'   Xmat or Bmat, except 'eps' which is a number. A zero precision will imply
#'   that the term is fixed.
#' @param init An initial values of the regression coefficients for the
#'   sampler...
#' @param burn The number of samples to be burned before actual sampling.
#' @param size The number of samples to be obtained from the samples.
#' 
#' @return A list of length two, posterior samples and means are stored in
#'   'samples' and 'means' elements respectively. Within each elements, 'coef'
#'   are the regression coefficients and 'prec' the precisions. The coef samples
#'   are organised as a list of arrays (dim_spl * length(level) * size) and a
#'   list of matrices (ncol(Xmat) * size) for effects. The prec samples are
#'   arrays for splines and effects terms, and vector for epsilon. The samples
#'   are always populated along the last dimension of the array/matrix/vector.
#'
bayes_ridge_semi_v4 <- function(y, Bmat, Xmat,
                                prior = NULL, prec = NULL, init = NULL,
                                burn = 0, size = 1000, debug = TRUE) {

  ## assumptions on y, Bmat, Xmat and grp. grp$sub are 'sticked' together,
  ## i.e. rows belonging to the same subjects are sticked together.

  check_Bmat(Bmat)
  check_Xmat(Xmat)
  which_sub <- purrr::map_lgl(Bmat, ~attr(.x, 'is_sub'))
  ## check_prec_v4(prec, Bmat, Xmat)

  ## put subject curves at the front
  Bmat <- Bmat[order(which_sub, decreasing = TRUE)]
  ## check if each element has a name
  stopifnot(length(setdiff(names(Bmat), '')) == length(Bmat))
  names(Bmat[[1]]) <- NULL # get rid of any potential names in the subject spline

  if (is.null(prec)) {
    prec_ls <- NULL
    stopifnot(!is.null(prior))
    stopifnot(c('spl', 'eff', 'eps') %in% names(prior),
              names(Bmat) %in% names(prior$spl),
              names(Xmat) %in% names(prior$eff))
    prior_ls <- list(spl = prior$spl[names(Bmat)],
                     eff = prior$eff[names(Xmat)],
                     eps = prior$eps)
    prior <- prec <- NULL
  } else {
    message('Emperical Bayes. Prior ignored.')
    prior_ls <- NULL
    stopifnot(c('spl', 'eff', 'eps') %in% names(prec),
              names(Bmat) %in% names(prec$spl),
              names(Xmat) %in% names(prec$eff))
    prec_ls <- list(spl = prec$spl[names(Bmat)],
                    eff = prec$eff[names(Xmat)],
                    eps = prec$eps)
    prior <- prec <- NULL
  }

  para <- get_para_v4(Bmat, Xmat)
  
  ## para: n_subs, subs_in_pop, pop_of_subs
  coef_samples <- get_coef_container_v4(para, size)
  prec_samples <- get_prec_container_v4(para, size)
  ## kcoef <- initialise_v4(Bmat, Xmat, y)
  kprec <- initialise_prec_v4(Bmat, Xmat, y)

  ## some pre-calculation
  xBs <- calc_xBs_v4(Bmat[[1]])
  
  ## ## make sure the names of precision/prior are in the same order as kcoef
  ## stopifnot(names(prior_ls$spl) == names(kcoef$spl),
  ##           names(prior_ls$eff) == names(kcoef$eff))    
  ## stopifnot(names(prec_ls$spl) == names(kcoef$spl),
  ##           names(prec_ls$eff) == names(kcoef$eff))  

  ## make sure the names of precision/prior are in the same order as kprec
  stopifnot(names(prior_ls$spl) == names(kprec$spl),
            names(prior_ls$eff) == names(kprec$eff))    
  stopifnot(names(prec_ls$spl) == names(kprec$spl),
            names(prec_ls$eff) == names(kprec$eff))  
  
  for (k in seq.int(-burn + 1, size)) {
    kcoef <- update_coefs_v4(y, Bmat, Xmat, xBs, kprec, debug)
    kprec <- update_precs_v4(kcoef, y, prior_ls, prec_ls)
    
    ## print progress
    if (isTRUE(k %% floor(size / 5) == 0)) {
      message(k, " samples generated.")
    }

    if (k > 0) {
      for (l in names(coef_samples$spl)) {
        stopifnot(isTRUE(all(dimnames(coef_samples$spline[[l]])[[2]] ==
                               colnames(kcoef$spl[[l]]))))
        coef_samples$spline[[l]][, , k] <- kcoef$spl[[l]]
        prec_samples$spline[[l]][, , k] <- kprec$spl[[l]]
      }
      for (l in names(coef_samples$eff)) {
        coef_samples$effect[[l]][, k] <- kcoef$eff[[l]]
        prec_samples$effect[[l]][, , k] <- kprec$eff[[l]]
      }
      prec_samples$eps[k] <- kprec$eps
    }
  }
  
  list(samples = list(coef = coef_samples,
                      prec = prec_samples),
       means = list(coef = recurse(coef_samples, pmean_v4),
                    prec = recurse(prec_samples, pmean_v4)))
}
