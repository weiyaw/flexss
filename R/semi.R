tcrossprod <- Matrix::tcrossprod
crossprod <- Matrix::crossprod

## inverse of a symmetrical, positive definite matrix
pdinv <- function(x) {
  stopifnot(Matrix::isSymmetric(x))
  x <- as(Matrix::forceSymmetric(x), 'dpoMatrix')
  Matrix::solve(x)
}

## product of the quadratic, x %*% A %*% t(x)
tquadprod <- function(x, A) {
  res <- tcrossprod(x %*% A, x)
  Matrix::forceSymmetric(res)
}

## the rows of Kmat must be independent
check_Kmat <- function(Kmat) {
  if (abs(Matrix::det(Matrix::tcrossprod(Kmat))) < 1e-10) {
    warning('The rows of Kmat may not be independent.')
  }
}

check_prec_beta <- function(prec_beta, dim_beta) {
  if (!is.null(prec_beta)) {
    message('prec$beta specified.')
    stopifnot(is.matrix(prec_beta),
              dim(prec_beta) == dim_beta,
              is_precision(prec_beta))
  }
}
check_prec <- function(prec, para) {
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



## is this a precision matrix (in the context of this model)? The precision can
## be semi-positive definite.
is_precision <- function(x) {
  all(isSymmetric(x), eigen(x, TRUE)$values >= 0)
}

init_theta <- function(theta, pop_names, dim_theta, n_pops, pls) {

  ## use pls if initial values not supplied
  if (is.null(theta)) {
    theta <- matrix(pls, dim_theta, n_pops)
  } else {
    message("Initial theta supplied.")
    if (!all(is.matrix(theta), dim(theta) == c(dim_theta, n_pops))) {
      stop("Invalid theta dimension. Must be a matrix")
    }
  }

  ## make sure the column names are correct
  if (is.null(colnames(theta))) {
    colnames(theta) <- pop_names
    theta
  } else {
    theta[, pop_names]
  }
}

init_delta <- function(delta, sub_names, dim_delta, n_subs) {

  ## use Gaussian rv if initial values not supplied
  if (is.null(delta)) {
    gauss_rv <- rnorm(dim_delta * n_subs) * 0.01
    delta <- matrix(gauss_rv, dim_delta, n_subs)
  } else {
    message("Initial delta supplied.")
    if (!all(is.matrix(delta), dim(delta) == c(dim_delta, n_subs))) {
      stop("Invalid delta dimension. Must be a matrix.")
    }
  }

  ## make sure the column names are correct
  if (is.null(colnames(delta))) {
    colnames(delta) <- sub_names
    delta
  } else {
    delta[, sub_names]
  }
}

## beta can be NULL, named (ordinary or column) vector
init_beta <- function(beta, beta_names, dim_beta) {
  
  ## use Gaussian rv if initial values not supplied
  if (is.null(dim_beta)) {
    return(NULL)
  } else if (is.null(beta)) {
    beta <- rnorm(dim_beta) * 0.01
  } else {
    message("Initial beta supplied.")
    if (!any(is.vector(beta, mode = 'numeric'),
             all(is.matrix(beta), dim(beta) == c(dim_beta, 1)))) {
      stop("Invalid beta dimension. Must be a vector/column vector.")
    }
  }
  ## make sure beta is a column vector
  beta <- as.matrix(beta)

  ## make sure the row names are correct
  if (is.null(rownames(beta))) {
    rownames(beta) <- beta_names
    beta
  } else {
    beta[beta_names, ]
  }
}

## datals: list(Bmat_pop: matrix, Bmat_sub: list of matrix, Xmat: matrix)
## theta, beta: col matrix
## delta: matrix
## return: col matrix
init_yhat <- function(datals, theta, delta, beta, debug = FALSE) {
  
  if (debug) {
  }
  f_init <- calc_f(datals$Bmat_pop, theta)
  g_init <- purrr::imap(datals$Bmat_sub, ~as.vector(.x %*% delta[, .y]))
  Xb_init <- calc_Xb(datals$Xmat, beta)

  f_init + unlist(g_init, use.names = FALSE) + Xb_init
}

## initialise theta with a penalised LS estimate, and delta with rnorms with small sd
initialise_with_pls_v3 <- function(init, para, pls, datals, debug = FALSE) {

  ## init: NULL or list(theta, delta, beta); any of the element can be NULL

  ## theta_names <- names(para$subs_in_pop)
  ## delta_names <- names(para$pop_of_subs)

  ## theta <- init_theta(init$theta, theta_names, para$dim_theta, para$n_pops, pls)
  theta <- as.matrix(rnorm(para$dim_theta)) # just for prototyping

  ## delta <- init_delta(init$delta, delta_names, para$dim_delta, para$n_subs)
  delta <- matrix(rnorm(para$dim_delta * para$n_subs),
                  para$dim_delta, para$n_subs,
                  dimnames = list(NULL, names(datals$Bmat_sub))) # just for prototyping

  ## beta <- init_beta(init$beta, para$beta_names, para$dim_beta)
  beta <- if(is.null(para$dim_beta)) {
    NULL # just for prototyping
  } else {
    as.matrix(rnorm(para$dim_beta)) # just for prototyping
  }

  yhat <- init_yhat(datals, theta, delta, beta, debug)

  ## the elements should have the same format as the output from update_coefs
  list(theta = theta, delta = delta, beta = beta, yhat = yhat)
}


## Bmat: list of matrix
## return: list of matrix
calc_xBs <- function(Bmat_sub) {
  purrr::map(Bmat_sub, Matrix::crossprod)
}

## xBs: list of matrix
## kprec: list(sub: matrix, eps: numeric)
## return: list of matrix
calc_Lmat_inv <- function(xBs, kprec) {
  calc_Lmat_i_inv <- function(xBs_i) {
    Lmat_i <- xBs_i + kprec$delta / kprec$eps
    pdinv(Lmat_i) # L_i > 0
  }
  purrr::map(xBs, calc_Lmat_i_inv)
}

## Bmat_sub: list of matrix (with 'idx' attr)
## Lmat: list of matrix
## return: matrix
calc_Phi <- function(Bmat_sub, Lmat_inv) {
  calc_Phi_i <- function(Bmat_sub_i, Lmat_i_inv) {
    dim_Phi <- NROW(Bmat_sub_i)
    diag(dim_Phi) - tquadprod(Bmat_sub_i,  Lmat_i_inv)
  }
  ## Phi <- purrr::map2(Bmat_sub, Lmat, calc_Phi_i)
  ## attr(Phi, 'idx') <- attr(Bmat_sub, 'idx')
  Phi_ls <- purrr::map2(Bmat_sub, Lmat_inv, calc_Phi_i)
  ## Phi_ls
  Matrix::.bdiag(Phi_ls)
}

## Phi: matrix
## Bmat_pop: matrix
## return: matrix (must match the entries in Bmat_pop)
calc_PhixB <- function(Phi, Bmat_pop) {
  ## Bmat_pop_ls <- split.data.frame(Bmat_pop, grp_sub)
  ## PhixB_ls <- purrr::map2(Phi, Bmat_pop_ls, `%*%`) # Phi_i is symmetrical
  ## do.call(rbind, PhixB_ls)
  Phi %*% Bmat_pop
}

## Bmat_pop: matrix
## PhixB: matrix
## kprec: list(theta: matrix, eps: numeric)
## return: matrix
calc_Nmat_inv <- function(Bmat_pop, PhixB, kprec) {
  Nmat <- Matrix::crossprod(Bmat_pop, PhixB) + kprec$theta / kprec$eps
  pdinv(Nmat) # N > 0
}

## Phi: matrix
## PhixB: matrix
## return: matrix
calc_Psy <- function(Phi, PhixB, Nmat_inv) {
  Phi - tquadprod(PhixB, Nmat_inv)
  ## calc_Psy_i <- function(Phi_i, PhixBs_i, Nmat_i_inv) {
  ##   Phi_i - tcrossprod(PhixBs_i %*% Nmat_i_inv, PhixBs_i) # Phi_i == t(Phi_i)
  ## }
  ## purrr::pmap(list(Phi, PhixBs, Nmat_inv), calc_Psy_i)
}

## Psy: matrix
## Bmat_pop: matrix
## return: matrix
calc_PsyxX <- function(Psy, Xmat) {
  if (is.null(Xmat)) {
    NULL
  } else {
    Psy %*% Xmat
  }
}


## Xmat: matrix
## PsyxX: matrix
## kprec: list(beta: matrix, eps: numeric)
## return: matrix
calc_Qmat_inv <- function(Xmat, PsyxX, kprec) {
  if (is.null(Xmat) || is.null(PsyxX)) {
    NULL
  } else {
    Qmat <- crossprod(Xmat, PsyxX) + kprec$beta / kprec$eps
    pdinv(Qmat)
  }
}

## y: col matrix
## PsyxX: matrix
## return: row matrix
calc_Rmat <- function(y, PsyxX) {
  if (is.null(PsyxX)) {
    NULL
  } else {
    crossprod(y, PsyxX)
  }
}

## Qmat_inv: matrix
## Rmat: matrix
## kprec: list(eps: numeric)
## return: col matrix
update_beta <- function(Qmat_inv, Rmat, kprec) {
  if (is.null(Qmat_inv) || is.null(Rmat)) {
    NULL
  } else {
    mu_beta <- as.vector(tcrossprod(Qmat_inv, Rmat))
    Sig_beta <- as.matrix(Qmat_inv / kprec$eps)
    t(mvtnorm::rmvnorm(1, mu_beta, Sig_beta))
  }
}

## Xmat: matrix
## beta: col matrix
## return: col matrix
calc_Xb <- function(Xmat, beta) {
  if (is.null(Xmat) || is.null(beta)) {
    NULL
  } else {
    Xmat %*% beta
  }
}

## y: col matrix
## Xb: col matrix
## PhixB: matrix
## return: row matrix
calc_Pmat <- function(y, Xb, PhixB) {
  if (is.null(Xb)) {
    crossprod(y, PhixB)
  } else {
    crossprod(y - Xb, PhixB)
  }
}

## Nmat_inv: matrix
## Pmat: row matrix
## kprec: list(eps: numeric)
## return: a list of length n_pop
update_theta <- function(Nmat_inv, Pmat, kprec) {
  mu_theta <- as.vector(tcrossprod(Nmat_inv, Pmat))
  Sig_theta <- as.matrix(Nmat_inv / kprec$eps)
  t(mvtnorm::rmvnorm(1, mu_theta, Sig_theta))
}

## Bmat_pop: matrix
## theta: col matrix
## return: col matrix
calc_f <- function(Bmat_pop, theta) {
  Bmat_pop %*% theta
}

## y: col matrix
## f: col matrix
## Xb: col matrix
## Bmat_sub: list of matrix
## grp_sub: factor
## return: list of row matrix
calc_Mmat <- function(y, f, Xb, Bmat_sub, grp_sub) {
  resids <- split(y - f - Xb, grp_sub)
  stopifnot(names(resids) == names(Bmat_sub))
  purrr::map2(resids, Bmat_sub, crossprod)
}

## Lmat_inv: list of matrix
## Mmat: list of matrix
## kprec: list(eps: numeric)
## return: list of col matrix
update_delta <- function(Lmat_inv, Mmat, kprec) {
  update_delta_i <- function(Lmat_inv_i, Mmat_i) {
    mu_delta <- as.vector(tcrossprod(Lmat_inv_i, Mmat_i))
    Sig_delta <- as.matrix(Lmat_inv_i / kprec$eps)
    t(mvtnorm::rmvnorm(1, mu_delta, Sig_delta))
  }
  purrr::map2(Lmat_inv, Mmat, update_delta_i)
}

## Bmat_sub: list of matrix
## delta: list of col matrix
## return: col matrix
calc_g <- function(Bmat_sub, delta) {
  g_ls <- purrr::map2(Bmat_sub, delta, ~as.vector(.x %*% .y))
  as.matrix(unlist(g_ls, use.names = FALSE))
}



## calc_xX <- function(Xmat) {
##   if (is.null(Xmat)) {
##     NULL
##   } else {
##     purrr::map(Xmat, crossprod)
##   }
## }

## ## return: list of length n_sub
## calc_Lmat <- function(xX, kprec, n_subs) {
##   if (is.null(xX)) {
##     NULL
##   } else if (all(kprec$beta == 0)) {
##     xX
##   } else {
##     purrr::map(xX, ~.x + kprec$beta / (n_subs * kprec$eps))
##   }
## }

## ## return: list of length n_sub
## calc_Mmat <- function(y, f, g, Xmat) {
##   if (is.null(Xmat)) {
##     NULL
##   } else {
##     purrr::pmap(list(y, f, g, Xmat), ~crossprod(..1 - ..2 - ..3, ..4))
##   }
## }

## ## return: list of length n_sub
## ## if Xmat is NULL (i.e. no beta term), y must be provided
## calc_Phi <- function(Xmat, Lmat, y = NULL) {
##   if (is.null(Xmat)) {
##     stopifnot(!is.null(y))
##     purrr::map(y, ~diag(length(.x)))
##   } else {
##     calc_Phi_i <- function(Xmat_i, Lmat_i) {
##       dim_Phi <- NROW(Xmat_i)
##       Lmat_i_inv <- chol2inv(chol(Lmat_i)) # L_i > 0
##       diag(dim_Phi) - tcrossprod(Xmat_i %*% Lmat_i_inv, Xmat_i)
##     }
##     purrr::map2(Xmat, Lmat, calc_Phi_i)
##   }
## }

## ## return: list of length n_sub
## calc_PhixBs <- function(Phi, Bmat_sub) {
##   purrr::map2(Phi, Bmat_sub, `%*%`) # Phi_i is symmetrical
## }

## ## return: list of length n_sub
## calc_Nmat_inv <- function(Bmat_sub, PhixBs, kprec) {
##   calc_Nmat_i_inv <- function(Bmat_sub_i, PhixBs_i) {
##     Nmat_i <- crossprod(Bmat_sub_i, PhixBs_i) + kprec$delta / kprec$eps
##     chol2inv(chol(Nmat_i)) # N_i > 0
##   }
##   purrr::map2(Bmat_sub, PhixBs, calc_Nmat_i_inv)
## }

## ## return: list of length n_sub
## calc_Pmat <- function(y, f, PhixBs) {
##   purrr::pmap(list(y, f, PhixBs), ~crossprod(..1 - ..2, ..3))
## }

## ## return: list of length n_sub
## calc_Psy <- function(Phi, PhixBs, Nmat_inv) {
##   calc_Psy_i <- function(Phi_i, PhixBs_i, Nmat_i_inv) {
##     Phi_i - tcrossprod(PhixBs_i %*% Nmat_i_inv, PhixBs_i) # Phi_i == t(Phi_i)
##   }
##   purrr::pmap(list(Phi, PhixBs, Nmat_inv), calc_Psy_i)
## }

## ## done
## ## return: list of length n_sub
## calc_PsyxB <- function(Psy, Bmat_pop) {
##   purrr::map2(Psy, Bmat_pop, `%*%`) # Psy_i is symmetrical
## }


## ## Bmat_pop, Psy : named list of length = n_sub
## ## subs_in_pop: named list of length = n_pop, each element a vector of sub names
## ## of that particular pop
## ## return: list of length n_pop
## calc_Qmat_inv <- function(Bmat_pop, PsyxB, subs_in_pop, kprec) {
##   ## require: magrittr
##   calc_Qmat_l_inv <- function(subs_l) {
##     subs_l <- as.character(subs_l)
##     Bmat_pop_l <- Bmat_pop[subs_l]
##     PsyxB_l <- PsyxB[subs_l]
##     Qmat_l <- purrr::map2(Bmat_pop_l, PsyxB_l, crossprod) %>%
##       purrr::reduce(`+`, .init = kprec$theta / kprec$eps) # this step gives Q_l
##     chol2inv(chol(Qmat_l))                             # Q_l > 0
##   }
##   purrr::map(subs_in_pop, calc_Qmat_l_inv)
## }

## ## return: list of length n_pop
## calc_Rmat <- function(y, PsyxB, subs_in_pop) {
##   calc_Rmat_l <- function(subs_l) {
##     subs_l <- as.character(subs_l)
##     y_l <- y[subs_l]
##     PsyxB_l <- PsyxB[subs_l]
##     purrr::map2(y_l, PsyxB_l, crossprod) %>%
##       purrr::reduce(`+`)
##   }
##   purrr::map(subs_in_pop, calc_Rmat_l)
## }

## ## pop_of_subs: named vector/list of length = n_sub of the pop of every sub
## ## return: list of length n_sub
## calc_f <- function(Bmat_pop, theta, pop_of_subs, debug = FALSE) {
##   if (debug) stopifnot(names(Bmat_pop) == names(pop_of_subs))
##   purrr::map2(Bmat_pop, pop_of_subs, ~.x %*% theta[[.y]])
## }

## ## return: list of length n_sub
## calc_g <- function(Bmat_sub, delta, debug = FALSE) {
##   if (debug) stopifnot(names(Bmat_sub) == names(delta))
##   purrr::map2(Bmat_sub, delta, `%*%`)
## }

## ## return: a list of length n_pop
## update_theta <- function(Bmat_pop, PsyxB, y, subs_in_pop, kprec) {
##   Qmat_inv <- calc_Qmat_inv(Bmat_pop, PsyxB, subs_in_pop, kprec)
##   Rmat <- calc_Rmat(y, PsyxB, subs_in_pop) 
##   mu_theta <- purrr::map2(Qmat_inv, Rmat, tcrossprod)
##   purrr::map2(mu_theta, Qmat_inv, ~t(mvtnorm::rmvnorm(1, .x, .y / kprec$eps)))
## }

## ## return: a list of length n_sub
## update_delta <- function(Nmat_inv, PhixBs, y, f, kprec) {
##   Pmat <- calc_Pmat(y, f, PhixBs)                # y, f
##   mu_delta <- purrr::map2(Nmat_inv, Pmat, tcrossprod)
##   purrr::map2(mu_delta, Nmat_inv, ~t(mvtnorm::rmvnorm(1, .x, .y / kprec$eps)))
## }

## ## return: a column vector
## update_beta <- function(Lmat, Mmat, y, f, g, kprec) {
##   if (is.null(Lmat) || is.null(Mmat)) {
##     NULL
##   } else {
##     Lmat_sum <- purrr::reduce(Lmat, `+`)
##     Lmat_sum_inv <- chol2inv(chol(Lmat_sum))
##     Mmat_sum <- purrr::reduce(Mmat, `+`)
##     mu_beta <- tcrossprod(Lmat_sum_inv, Mmat_sum)
##     Sig_beta <- Lmat_sum_inv / kprec$eps
##     t(mvtnorm::rmvnorm(1, mu_beta, Sig_beta))
##   }
## }

## ## beta: a vector / column vector of length = NCOL(Xmat)
## ## the rest of args: list of column vectors/matrix, length = n_subs
## ## return a column vector (not numeric vector)
## update_yhat <- function(f, g, Xmat, beta) {
##   ## stopifnot(names(f) == names(g))
##   ## if (is.null(Xmat)) {
##   ##   unlist(f) + unlist(g)
##   ## } else {
##   ##   stopifnot(names(f) == names(Xmat))
##   ##   tall_Xmat <- do.call(rbind, Xmat)
##   ##   unlist(f) + unlist(g) + tall_Xmat %*% beta
##   ## }

##   ## safer to produce a list, or risk messing up the subject names
##   stopifnot(names(f) == names(g))
##   if (is.null(Xmat)) {
##     purrr::map2(f, g, `+`)
##   } else {
##     stopifnot(names(f) == names(Xmat))
##     purrr::pmap(list(f, g, Xmat), ~..1 + ..2 + ..3 %*% beta)
##   }
## }


## datals: list(Bmat_pop: matrix, Bmat_sub: list of matrix, Xmat: matrix, y: col matrix)
## xBs: list of matrix
## para: list(grp_sub: factor)
## kprec: list(theta: matrix, delta: matrix, beta: matrix, eps: numeric)
update_coefs <- function(datals, xBs, para, kprec, debug = FALSE) {

  if (debug) {
    ## stopifnot(names(datals$Bmat_pop) == names(datals$Bmat_sub),
    ##           names(datals$Bmat_pop) == names(datals$Xmat),
    ##           names(datals$Bmat_pop) == names(datals$y),
    ##           names(datals$Bmat_pop) == names(xX),
    ##           names(datals$Bmat_pop) == names(para$pop_of_subs))
  }
  
  ## some pre-calculation
  Lmat_inv <- calc_Lmat_inv(xBs, kprec[c('delta', 'eps')])
  Phi <- calc_Phi(datals$Bmat_sub, Lmat_inv)
  PhixB <- calc_PhixB(Phi, datals$Bmat_pop) # Phi_i is symmetrical
  Nmat_inv <- calc_Nmat_inv(datals$Bmat_pop, PhixB, kprec[c('theta', 'eps')])
  Psy <- calc_Psy(Phi, PhixB, Nmat_inv)
  PsyxX <- calc_PsyxX(Psy, datals$Xmat) # Psy_i is symmetrical

  ## update beta
  Qmat_inv <- calc_Qmat_inv(datals$Xmat, PsyxX, kprec[c('beta', 'eps')])
  Rmat <- calc_Rmat(datals$y, PsyxX)
  beta <- update_beta(Qmat_inv, Rmat, kprec['eps'])
  Xb <- calc_Xb(datals$Xmat, beta)
    
  ## update theta
  Pmat <- calc_Pmat(datals$y, Xb, PhixB)
  theta <- update_theta(Nmat_inv, Pmat, kprec['eps'])
  f <- calc_f(datals$Bmat_pop, theta)
  
  ## update delta
  Mmat <- calc_Mmat(datals$y, f, Xb, datals$Bmat_sub, para$grp_sub)
  delta <- update_delta(Lmat_inv, Mmat, kprec)
  g <- calc_g(datals$Bmat_sub, delta)
  
  ## delta: matrix, each column represents one /subject
  ## theta, beta: column vectors
  ## yhat: a list of column vectors, length = n_subs
  list(theta = theta,
       delta = do.call(cbind, delta) %>% `colnames<-`(names(delta)),
       beta = beta,
       yhat = f + g + Xb)
}

get_coef_container <- function(para, size) {
  ## theta_array <- array(NA, c(para$dim_theta, para$n_pops, size),
  ##                      dimnames = list(NULL, names(para$subs_in_pop)))
  theta_array <- matrix(NA, para$dim_theta, size) # just for prototyping
  delta_array <- array(NA, c(para$dim_delta, para$n_subs, size),
                       dimnames = list(NULL, unique(para$grp_sub)))
  if (is.null(para$dim_beta)) {
    beta_matrix <- NULL
  } else {
    beta_matrix <- matrix(NA, para$dim_beta, size,
                          dimnames = list(para$beta_names))
  }
  list(theta = theta_array,
       delta = delta_array,
       beta = beta_matrix)
}

get_prec_container <- function(para, size, init_prec) {
  theta_array <- array(NA, c(para$dim_theta, para$dim_theta, size))
  delta_array <- array(NA, c(para$dim_delta, para$dim_delta, size))
  eps_vector <- rep(NA, size)

  if (is.null(para$dim_beta)) {
    beta_array <- NULL
  } else {
    beta_array <- array(NA, c(para$dim_beta, para$dim_beta, size))
  }

  list(theta = theta_array,
       delta = delta_array,
       beta = beta_array,
       eps = eps_vector)
}


get_dims <- function(Bmat, Xmat) {
  list(theta = NCOL(Bmat$pop),
       delta = NCOL(Bmat$sub),
       beta = NCOL(Xmat))
}

get_para <- function(grp, Bmat, Xmat, Kmat, dim_block, ranef, prec_beta) {
  para <- list(n_pops = length(unique(grp$pop)),
               n_subs = length(unique(grp$sub)),
               dim_theta = NCOL(Bmat$pop),
               dim_delta = NCOL(Bmat$sub),
               dim_beta = NCOL(Xmat),
               Kmat = Kmat,
               dim_block = dim_block,
               ranef = ranef)
  para$grp_sub <- grp$sub
  ## para$subs_in_pop <- tapply(grp$sub, grp$pop,
  ##                            function(x) as.character(unique(x)),
  ##                            simplify = FALSE)
  ## para$pop_of_subs <- tapply(grp$pop, grp$sub,
  ##                            function(x) as.character(unique(x)),
  ##                            simplify = FALSE)
  
  if (is.null(Xmat)) {
    ## for ranef to be NULL if no random/fixed effects
    para[c('dim_beta', 'ranef', 'beta_names')] <- NULL
  } else {
    if (!is.null(prec_beta)) {
      message("ranef overriden by prec$beta") 
      para$ranef <- which(diag(prec_beta) > 0)
    }
    if (is.null(colnames(Xmat))) {
      para$beta_names <- paste0('beta', seq_len(NCOL(Xmat)))
    } else {
      para$beta_names <- colnames(Xmat)
    }
  }
  para
}

## x: a vector or matrix, to be squared and sum
update_with_gamma <- function(x, a, b) {
  stopifnot(!is.null(x), !is.null(a), !is.null(b))
  shape <- 0.5 * length(x) + a
  rate <- 0.5 * sum(x^2) + b
  rgamma(1, shape = shape, rate = rate)
}

## x: a matrix of col vectors x, to be tcrossprod
update_with_wishart <- function(x, v, lambda) {
  stopifnot(!is.null(x), !is.null(v), !is.null(lambda))
  df <- v + NCOL(x)
  scale <- lambda + tcrossprod(x)
  inv_scale <- pdinv(scale)
  rWishart(1, df = df, Sigma = inv_scale)[, , 1]
}

## theta: matrix, each column the coef of a population
## Kmat: the penalty matrix (row vectors must be independent)
## Kmat %*% theta are penalised (i.e. are random)
## return precision matrix for theta
update_prec_theta <- function(theta, Kmat, prior) {
  x <- Kmat %*% theta
  update_with_gamma(x, prior$a, prior$b) * Matrix::crossprod(Kmat)
}

## delta: matrix, each column the coef of a subject
## dim_block: the dimension corresponding to the block cov matrix
## return precision matrix for delta
update_prec_delta <- function(delta, dim_block, prior) {

  block <- NA
  iid <- NA
  
  if (dim_block > NROW(delta) || dim_block < 0) {
    stop('dim_block out of bound.')
  }
  
  ## block precision term
  if (dim_block > 0) {
    delta_block <- delta[seq(1, dim_block), , drop = FALSE]
    block <- update_with_wishart(delta_block, prior$v, prior$lambda)
  }
  ## iid precision term
  if (dim_block < NROW(delta)) {
    delta_iid <- delta[seq(dim_block + 1, NROW(delta)), , drop = FALSE]
    iid <- update_with_gamma(delta_iid, prior$a, prior$b)
  }

  diag(iid, NROW(delta)) %>%
    `[<-`(seq_len(dim_block), seq_len(dim_block), block)
}

## yhat, y: list of col vectors/vectors with the same names
## return precision of the residual
update_prec_eps <- function(yhat, y, prior) {
  stopifnot(names(y) == names(yhat))
  resids <- y - yhat
  ## resids <- unlist(y) - unlist(yhat)
  update_with_gamma(resids, prior$a, prior$b)
}

## beta: col vector/vector
## ranef: a vector of indices that corresponds to the random effects
## if ranef is NULL, return a zero matrix. (i.e. everything is fixed effect)
update_prec_beta <- function(beta, ranef, prior) {
  if (is.null(beta)) {
    NULL
  } else if (is.null(ranef)) {
    diag(0, length(beta))
  } else {
    stopifnot(NCOL(beta) == 1)
    x <- beta[ranef]
    prec <- update_with_gamma(x, prior$a, prior$b)
    rep(0, times = NROW(beta)) %>%
      `[<-`(ranef, prec) %>%
      diag()
  }
}

update_precs <- function(kcoef, y, para, prior_ls, init_prec = NULL) {
  if (!is.null(init_prec)) {
    if (is.null(para$dim_beta)) {
      init_prec$beta <- NULL
    }
    init_prec
  } else {
    list(theta = update_prec_theta(kcoef$theta, para$Kmat, prior_ls$theta),
         delta = update_prec_delta(kcoef$delta, para$dim_block, prior_ls$delta),
         eps = update_prec_eps(kcoef$yhat, y, prior_ls$eps),
         beta = update_prec_beta(kcoef$beta, para$ranef, prior_ls$beta))
  }
}


## Bmat: list(pop: matrix, sub: matrix)
## Xmat: matrix
## y: col matrix
## para: list(grp_sub: factor)
## split data into chucks according to the subject index arguments are same as
## bayes_ridge_semi().
split_data_subjects <- function(Bmat, Xmat, y, para, debug) {

  stopifnot(is.factor(para$grp_sub)) # must be a factor, or risk mixing up row index

  Bmat_sub <- split.data.frame(Bmat$sub, para$grp_sub)

  datals <- list(Bmat_pop = Matrix::Matrix(Bmat$pop, sparse = TRUE),
                 Bmat_sub = purrr::map(Bmat_sub, Matrix::Matrix),
                 y = Matrix::Matrix(y))

  if (!is.null(Xmat)) {
    datals$Xmat <- Matrix::Matrix(Xmat, sparse = TRUE)
  }

  
  if (debug) {
    ## stopifnot(names(datals$Bmat_pop) == names(datals$Bmat_sub),
    ##           names(datals$Bmat_pop) == names(datals$Xmat),
    ##           names(datals$Bmat_pop) == names(datals$y))
  }
  datals
}

## calculate posterior mean from an array or vector of samples. The samples are
## populated along the last dimension, e.g. the columns of a matrix are the
## samples; the depth of a 3D array are the samples.
pmean <- function(samples) {
  if (is.array(samples)) {
    rowMeans(samples, dims = length(dim(samples)) - 1)
  } else if (is.vector(samples)) {
    mean(samples)
  } else {
    stop("Invalid samples structure.")
  }
}

## This is a Gibbs sampler v3 for longitudinal Bayesian semiparametric ridge.
## Update two block of parameters: variance and coefs

## Requirements: y, grp, Bmat (df), Xmat (df), Kmat, dim_block, ranef

## Algorithms paremeters: burn, size
## Extras: init (list of matrices with 'pop' and 'sub'), prior (see 'check_prior')


#' Bayesian ridge for longitudinal semiparemetric models
#'
#' This is a Gibbs sampler v3 for fitting Bayesian longitudinal semiparametric
#' models. It assumes that there are multiple populations in the model and each
#' population is modelled by a 'population' curve. Within each population, they
#' are multiple subjects, and each subject are modelled by a 'subject'
#' curve. The subject curves are treated as deviations from their respective
#' population curves. On top of that, fixed or random effects can be added to
#' the model. This is useful when, for example, a particular treatment was
#' applied to some of the populations or subjects, and the user is interested in
#' the effect of the treatment.
#'
#' The population and subject curves are modelled as linear (in statistical
#' sence, not only stright lines). This includes polynomials and splines.
#'
#' The mathematical model is
#'
#' y_i = Bmat$pop_i %*% \theta_i + Bmat$sub_i %*% \delta_i + Xmat_i %*% \beta +
#' \epsilon_i
#'
#' where 'i' is the index of subjects, y_i is a vector of observation for the
#' i^{th} subjectm, Bmat$pop_i and Bmat$sub_i are the model/design matrices for
#' the population and subject curves (of the i^{th} subject) respectively, and
#' \epsilon_i are Gaussian errors. The \theta_i should be intepreted as the
#' regression coefficients of the population curve of the i^{th} subject, and
#' the subjects belonging to the same population will have the same \theta. This
#' implies that some of the \theta_i will be identical. However, \delta_i will
#' be different for each subject, while \beta will be the same for all subjects.
#'
#' The coefficients \theta, \delta, \beta and \epsilon are all Gaussian, but
#' their covariance structures are not necessarily diagonal. Consult the
#' manual/paper for more info.
#'
#' The sampler first updates the precision, then the coefficients. The joint
#' conditional posterior of all the coefficients is often highly correlated, but
#' fortunately a closed-form expression is available and is utilised here. This
#' implies that the convergence of this Gibbs sampler is quick and does not
#' require burning many samples in the beginning.
#' 
#' @param y A vector of the response vector.
#' @param grp A list of two numeric/factor/character vectors, grp$pop and
#'   grp$sub, specifying the population and subject to which each observation
#'   belongs.
#' @param Bmat A list of two design matrices for the splines, Bmat$pop and
#'   Bmat$sub for population- and subject-specific curves respectively.
#' @param Xmat A design matrix for the non-spline terms, i.e. fixed and random
#'   effects.
#' @param Kmat The penalty matrix applied on the population splines
#'   coefficients.
#' @param dim_block The number of subject spline coefficients to be treated as
#'   correlated. The first `dim_block` of the coefficients are correlated. Set
#'   this to 0 if all subject coefficients are uncorrelated, or `dim(Bmat$sub)`
#'   if everything is correlated.
#' @param ranef A numeric vector of indices of the Xmat columns where the
#'   corresponding coefficients should be treated as random effects. NULL if
#'   everything is fixed. This option will get overridden by the 'beta' term in
#'   `prec`.
#' @param burn The number of samples to be burned before actual sampling.
#' @param size The number of samples to be obtained from the samples.
#' @param init An initial values of the regression coefficients for the
#'   sampler...
#' @param prior unknown yet
#' @param prec An optional list of precisions with 'theta', 'delta', 'beta' and
#'   'eps'. If supplied, the precisions of the model are fixed. All the elements
#'   are square matrices, except 'eps' which is a number. The 'beta' term here
#'   will override `ranef` (i.e. all the 0 in diagonal will result in the
#'   corresponding beta being fixed).
#'
#' @return A list of length two, posterior samples and means are stored in
#'   'samples' and 'means' elements respectively. Within each elements, 'coef'
#'   are the regression coefficients and 'prec' the precisions (i.e. the inverse
#'   of variances).
#'
#'  The samples are organised as arrays, except coef$beta (matrix) and prec$eps
#'  (vector). The samples are always populated along the last dimension of the
#'  array, matrix or vector. For the 'coef' samples, the first dimension of the
#'  array/matrix is always the dimension of the coefficients, followed by the
#'  indices of populations/subjects (for 'theta' and 'delta' only). All the
#'  precisions except 'eps' are symmetrical square matrices.
#' 
bayes_ridge_semi <- function(y, grp, Bmat, Xmat,
                             Kmat = NULL, dim_block = NULL, ranef = NULL,
                             burn = 0, size = 1000, init = NULL, prior = NULL,
                             prec = NULL, debug = TRUE) {
  
  ## sort the data according to the subject index
  idx <- order(grp$sub, grp$pop) # grp$sub should be nested within grp$pop, so
                                 # the second argument should not matter
  y <- y[idx]
  grp <- purrr::map(grp, ~factor(.x[idx], unique(.x[idx])))
  Bmat <- purrr::map(Bmat, ~.x[idx, ])
  Xmat <- Xmat[idx, ]
  
  check_Kmat(Kmat)
  para <- get_para(grp, Bmat, Xmat, Kmat, dim_block, ranef, prec$beta)
  ## check_prec(prec, para)
  ## grp <- check_grp(grp)
  
  ## para: n_subs, subs_in_pop, pop_of_subs
  coef_samples <- get_coef_container(para, size)
  prec_samples <- get_prec_container(para, size, prec)

  ## get a list(Bmat_pop: matrix, Bmat_sub: list of matrix, y: col matrix, Xmat: matrix)
  datals <- split_data_subjects(Bmat, Xmat, y, para['grp_sub'], debug)

  pls <- get_pls(y, Bmat$pop, Kmat)
  kcoef <- initialise_with_pls_v3(init, para, pls, datals, debug)

  ## some pre-calculation
  xBs <- calc_xBs(datals$Bmat_sub)
  prior_ls <- prior
  
  for (k in seq.int(-burn + 1, size)) {
    kprec <- update_precs(kcoef, datals$y, para, prior_ls, prec)
    kcoef <- update_coefs(datals, xBs, para, kprec, debug)

    ## print progress
    if (k %% 1000 == 0) {
      message(k, " samples generated.")
    }

    if (k > 0) {
      coef_samples$theta[, k] <- kcoef$theta
      coef_samples$delta[, colnames(kcoef$delta), k] <- kcoef$delta
      coef_samples$beta[, k] <- kcoef$beta
      
      prec_samples$theta[, , k] <- as.matrix(kprec$theta)
      prec_samples$delta[, , k] <- as.matrix(kprec$delta)
      prec_samples$beta[, , k] <- as.matrix(kprec$beta)
      prec_samples$eps[k] <- kprec$eps
    }
  }
  
  list(samples = list(coef = coef_samples,
                      prec = prec_samples),
       means = list(coef = purrr::map(coef_samples, pmean),
                    prec = purrr::map(prec_samples, pmean)))
}
