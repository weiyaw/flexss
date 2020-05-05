
## the rows of Kmat must be independent
check_Kmat <- function(Kmat) {
  if (abs(det(tcrossprod(Kmat))) < 1e-10) {
    warning('The rows of Kmat may not be independent.')
  }
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
  if (is.null(beta)) {
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

## datals: list of lists
## theta, delta: matrix
## beta: column vector
## pop_of_subs: numeric vector/list of numbers
## return list of column vectors
init_yhat <- function(datals, theta, delta, beta, pop_of_subs, debug = FALSE) {

  if (debug) {
    stopifnot(names(datals$Bmat_pop) == names(pop_of_subs))
  }
  f_init <- purrr::map2(datals$Bmat_pop, pop_of_subs, ~.x %*% theta[, .y])
  g_init <- purrr::imap(datals$Bmat_sub, ~.x %*% delta[, .y])
  calc_yhat(f_init, g_init, datals$Xmat, beta)
}

## initialise theta with a penalised LS estimate, and delta with rnorms with small sd
initialise_with_pls_v3 <- function(init, para, pls, datals, debug = FALSE) {

  ## init: NULL or list(theta, delta, beta); any of the element can be NULL

  theta_names <- names(para$subs_in_pop)
  delta_names <- names(para$pop_of_subs)

  theta <- init_theta(init$theta, theta_names, para$dim_theta, para$n_pops, pls)
  delta <- init_delta(init$delta, delta_names, para$dim_delta, para$n_subs)
  beta <- init_beta(init$beta, para$beta_names, para$dim_beta)
  yhat <- init_yhat(datals, theta, delta, beta, para$pop_of_subs, debug)

  ## the elements should have the same format as the output from update_coefs
  list(theta = theta, delta = delta, beta = beta, yhat = yhat)
}


## return: list of length n_sub
calc_Lmat <- function(xX, kprec, n_subs) {
  purrr::map(xX, ~.x + kprec$beta / (n_subs * kprec$eps))
}

## return: list of length n_sub
calc_Mmat <- function(y, f, g, Xmat) {
  purrr::pmap(list(y, f, g, Xmat), ~crossprod(..1 - ..2 - ..3, ..4))
}

## return: list of length n_sub
calc_Phi <- function(Xmat, Lmat) {
  calc_Phi_i <- function(Xmat_i, Lmat_i) {
    dim_Phi <- NROW(Xmat_i)
    Lmat_i_inv <- chol2inv(chol(Lmat_i)) # L_i > 0
    diag(dim_Phi) - tcrossprod(Xmat_i %*% Lmat_i_inv, Xmat_i)
  }
  purrr::map2(Xmat, Lmat, calc_Phi_i)
}

## return: list of length n_sub
calc_PhixBs <- function(Phi, Bmat_sub) {
  purrr::map2(Phi, Bmat_sub, `%*%`) # Phi_i is symmetrical
}

## return: list of length n_sub
calc_Nmat_inv <- function(Bmat_sub, PhixBs, kprec) {
  calc_Nmat_i_inv <- function(Bmat_sub_i, PhixBs_i) {
    Nmat_i <- crossprod(Bmat_sub_i, PhixBs_i) + kprec$delta / kprec$eps
    chol2inv(chol(Nmat_i)) # N_i > 0
  }
  purrr::map2(Bmat_sub, PhixBs, calc_Nmat_i_inv)
}

## return: list of length n_sub
calc_Pmat <- function(y, f, PhixBs) {
  purrr::pmap(list(y, f, PhixBs), ~crossprod(..1 - ..2, ..3))
}

## return: list of length n_sub
calc_Psy <- function(Phi, PhixBs, Nmat_inv) {
 calc_Psy_i <- function(Phi_i, PhixBs_i, Nmat_i_inv) {
    Phi_i - tcrossprod(PhixBs_i %*% Nmat_i_inv, PhixBs_i) # Phi_i == t(Phi_i)
  }
  purrr::pmap(list(Phi, PhixBs, Nmat_inv), calc_Psy_i)
}

## done
## return: list of length n_sub
calc_PsyxB <- function(Psy, Bmat_pop) {
  purrr::map2(Psy, Bmat_pop, `%*%`) # Psy_i is symmetrical
}


## Bmat_pop, Psy : named list of length = n_sub
## subs_in_pop: named list of length = n_pop, each element a vector of sub names
## of that particular pop
## return: list of length n_pop
calc_Qmat_inv <- function(Bmat_pop, PsyxB, subs_in_pop, kprec) {
  ## require: magrittr
  calc_Qmat_l_inv <- function(subs_l) {
    subs_l <- as.character(subs_l)
    Bmat_pop_l <- Bmat_pop[subs_l]
    PsyxB_l <- PsyxB[subs_l]
    Qmat_l <- purrr::map2(Bmat_pop_l, PsyxB_l, crossprod) %>%
      purrr::reduce(`+`, .init = kprec$theta / kprec$eps) # this step gives Q_l
    chol2inv(chol(Qmat_l))                             # Q_l > 0
  }
  purrr::map(subs_in_pop, calc_Qmat_l_inv)
}

## return: list of length n_pop
calc_Rmat <- function(y, PsyxB, subs_in_pop) {
  calc_Rmat_l <- function(subs_l) {
    subs_l <- as.character(subs_l)
    y_l <- y[subs_l]
    PsyxB_l <- PsyxB[subs_l]
    purrr::map2(y_l, PsyxB_l, crossprod) %>%
      purrr::reduce(`+`)
  }
  purrr::map(subs_in_pop, calc_Rmat_l)
}

## pop_of_subs: named vector/list of length = n_sub of the pop of every sub
## return: list of length n_sub
calc_f <- function(Bmat_pop, theta, pop_of_subs, debug = FALSE) {
  if (debug) stopifnot(names(Bmat_pop) == names(pop_of_subs))
  purrr::map2(Bmat_pop, pop_of_subs, ~.x %*% theta[[.y]])
}

## return: list of length n_sub
calc_g <- function(Bmat_sub, delta, debug = FALSE) {
  if (debug) stopifnot(names(Bmat_sub) == names(delta))
  purrr::map2(Bmat_sub, delta, `%*%`)
}

## beta: a vector / column vector of length = NCOL(Xmat)
## the rest of args: list of column vectors/matrix, length = n_subs
calc_yhat <- function(f, g, Xmat, beta) {
  purrr::pmap(list(f, g, Xmat), ~..1 + ..2 + ..3 %*% beta)
}

## return: a list of length n_pop
update_theta <- function(Bmat_pop, PsyxB, y, subs_in_pop, kprec) {
  Qmat_inv <- calc_Qmat_inv(Bmat_pop, PsyxB, subs_in_pop, kprec)
  Rmat <- calc_Rmat(y, PsyxB, subs_in_pop) 
  mu_theta <- purrr::map2(Qmat_inv, Rmat, tcrossprod)
  purrr::map2(mu_theta, Qmat_inv, ~t(mvtnorm::rmvnorm(1, .x, .y / kprec$eps)))
}
  
## return: a list of length n_sub
update_delta <- function(Nmat_inv, PhixBs, y, f, kprec) {
  Pmat <- calc_Pmat(y, f, PhixBs)                # y, f
  mu_delta <- purrr::map2(Nmat_inv, Pmat, tcrossprod)
  purrr::map2(mu_delta, Nmat_inv, ~t(mvtnorm::rmvnorm(1, .x, .y / kprec$eps)))
}

## return: a column vector
update_beta <- function(Lmat, Mmat, y, f, g, kprec) {
  Lmat_sum <- purrr::reduce(Lmat, `+`)
  Lmat_sum_inv <- chol2inv(chol(Lmat_sum))
  Mmat_sum <- purrr::reduce(Mmat, `+`)
  mu_beta <- tcrossprod(Lmat_sum_inv, Mmat_sum)
  Sig_beta <- Lmat_sum_inv / kprec$eps
  t(mvtnorm::rmvnorm(1, mu_beta, Sig_beta))
}

## datals (list of lists): Bmat_pop, Bmat_sub, Xmat, y
## xX (list)
## para (list): n_subs, subs_in_pop, pop_of_subs
## kprec (list)
update_coefs <- function(datals, xX, para, kprec, debug = FALSE) {

  if (debug) {
    stopifnot(names(datals$Bmat_pop) == names(datals$Bmat_sub),
              names(datals$Bmat_pop) == names(datals$Xmat),
              names(datals$Bmat_pop) == names(datals$y),
              names(datals$Bmat_pop) == names(xX),
              names(datals$Bmat_pop) == names(para$pop_of_subs))
  }

  ## some pre-calculation
  Lmat <- calc_Lmat(xX, kprec, para$n_subs)
  Phi <- calc_Phi(datals$Xmat, Lmat)
  PhixBs <- calc_PhixBs(Phi, datals$Bmat_sub)
  Nmat_inv <- calc_Nmat_inv(datals$Bmat_sub, PhixBs, kprec)

  ## update theta
  Psy <- calc_Psy(Phi, PhixBs, Nmat_inv)
  PsyxB <- calc_PsyxB(Psy, datals$Bmat_pop)
  theta <- update_theta(datals$Bmat_pop, PsyxB, datals$y, para$subs_in_pop, kprec)

  ## update delta
  f <- calc_f(datals$Bmat_pop, theta, para$pop_of_subs, debug)
  delta <- update_delta(Nmat_inv, PhixBs, datals$y, f, kprec)
  
  ## update beta (a column vector)
  g <- calc_g(datals$Bmat_sub, delta, debug)
  Mmat <- calc_Mmat(datals$y, f, g, datals$Xmat)
  beta <- update_beta(Lmat, Mmat, datals$y, f, g, kprec)
  
  ## calc prediction
  yhat <- calc_yhat(f, g, datals$Xmat, beta)

  ## theta & delta are matrices, each column represents one population/subject
  ## beta: a column vector
  ## yhat: a list of column vectors, length = n_subs
  list(theta = do.call(cbind, theta) %>% `colnames<-`(names(theta)),
       delta = do.call(cbind, delta) %>% `colnames<-`(names(delta)),
       beta = beta,
       yhat = yhat)
}
  
get_empty_coef_samples <- function(para, size) {
  theta_array <- array(NA, c(para$dim_theta, para$n_pops, size),
                     dimnames = list(NULL, names(para$subs_in_pop)))
  delta_array <- array(NA, c(para$dim_delta, para$n_subs, size),
                     dimnames = list(NULL, names(para$pop_of_subs)))
  beta_matrix <- matrix(NA, para$dim_beta, size,
                        dimnames = list(para$beta_names))
  list(population = theta_array,
       subjects = delta_array,
       beta = beta_matrix)
}

get_empty_prec_samples <- function(para, size) {
  theta_array <- array(NA, c(para$dim_theta, para$dim_theta, size))
  delta_array <- array(NA, c(para$dim_delta, para$dim_delta, size))
  beta_array <- array(NA, c(para$dim_beta, para$dim_beta, size))
  eps_vector <- rep(NA, size)

  list(theta = theta_array,
       delta = delta_array,
       beta = beta_array,
       eps = eps_vector)
}


get_para <- function(grp, Bmat, Xmat, Kmat, dim_block, ranef) {
  para <- list(n_pops = length(unique(grp$pop)),
               n_subs = length(unique(grp$sub)),
               dim_theta = NCOL(Bmat$pop),
               dim_delta = NCOL(Bmat$sub),
               dim_beta = NCOL(Xmat),
               Kmat = Kmat,
               dim_block = dim_block,
               ranef = ranef)
  para$subs_in_pop <- tapply(grp$sub, grp$pop,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  para$pop_of_subs <- tapply(grp$pop, grp$sub,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  
  if (is.null(colnames(Xmat))) {
    para$beta_names <- paste0('beta', seq_len(NCOL(Xmat)))
  } else {
    para$beta_names <- colnames(Xmat)
  }
  para
}

## x: a vector or matrix, to be squared and sum
update_with_gamma <- function(x, a, b) {
  shape <- 0.5 * length(x) + a
  rate <- 0.5 * sum(x^2) + b
  rgamma(1, shape = shape, rate = rate)
}

## x: a matrix of col vectors x, to be tcrossprod
update_with_wishart <- function(x, v, lambda) {
  df <- v + NCOL(x)
  scale <- lambda + tcrossprod(x)
  inv_scale <- chol2inv(chol(scale))
  rWishart(1, df = df, Sigma = inv_scale)[, , 1]
}

## theta: matrix, each column the coef of a population
## Kmat: the penalty matrix (row vectors must be independent)
## Kmat %*% theta are penalised (i.e. are random)
## return precision matrix for theta
update_prec_theta <- function(theta, Kmat, prior) {
  x <- Kmat %*% theta
  update_with_gamma(x, prior$a, prior$b) * crossprod(Kmat)
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
  resid <- unlist(yhat) - unlist(y)
  update_with_gamma(resid, prior$a, prior$b)
}

## beta: col vector/vector
## ranef: a vector of indices that corresponds to the random effects
update_prec_beta <- function(beta, ranef, prior) {
  stopifnot(NCOL(beta) == 1)
  x <- beta[ranef]
  prec <- update_with_gamma(x, prior$a, prior$b)
  rep(0, times = NROW(beta)) %>%
    `[<-`(ranef, prec) %>%
    diag()
}


update_precs <- function(kcoef, y, para, prior_ls) {
  list(theta = update_prec_theta(kcoef$theta, para$Kmat, prior_ls$theta),
       delta = update_prec_delta(kcoef$delta, para$dim_block, prior_ls$delta),
       eps = update_prec_eps(kcoef$yhat, y, prior_ls$eps),
       beta = update_prec_beta(kcoef$beta, para$ranef, prior_ls$beta))
}


## split data into chucks according to the subject index
## arguments are same as bayes_ridge_semi()
split_data_subjects <- function(Bmat, Xmat, y, grp_sub, debug) {
  datals <- list(Bmat_pop = split.data.frame(Bmat$pop, grp_sub),
                 Bmat_sub = split.data.frame(Bmat$sub, grp_sub),
                 Xmat = split.data.frame(Xmat, grp_sub),
                 y = split.default(y, grp_sub))

  if (debug) {
    stopifnot(names(datals$Bmat_pop) == names(datals$Bmat_sub),
              names(datals$Bmat_pop) == names(datals$Xmat),
              names(datals$Bmat_pop) == names(datals$y))
  }
  datals
}

## f1 <- function(n_samples, resids, ig_a, ig_b ) {
##   ## update sigma^2_eps
##   shape_eps <- 0.5 * n_samples + ig_a
##   rate_eps <- 0.5 * crossprod(resids) + ig_b
##   rgamma(1, shape = shape_eps, rate = rate_eps)
## }

## ## update sigma^2_theta
## shape_pop <- 0.5 * n_pops * rank_K + ig_a$pop
## rate_pop <- 0.5 * sum((Kmat %*% coef$pop)^2) + ig_b$pop
## prec$pop <- rgamma(1, shape = shape_pop, rate = rate_pop)

## f2 <- function() {
##   coef_sub <- matrix(1:6, 2, 3)
##   n_subs <- ncol(coef_sub)
##   iw_v <- 3
##   iw_lambda <- diag(nrow(coef_sub))
##   ## update Sigma_dev1
##   df_sub1 <- iw_v + n_subs
##   scale_sub1 <- iw_lambda + tcrossprod(coef_sub)
##   inv_scale_sub1 <- chol2inv(chol(scale_sub1))
##   rWishart(1, df = df_sub1, Sigma = inv_scale_sub1)[, , 1]
## }

## ## update sigma^2_dev2
## shape_sub2 <- 0.5 * n_subs * (dim_delta - dim_sub1) + ig_a$sub2
## rate_sub2 <- 0.5 * sum(coef$sub[-(1:dim_sub1), ]^2) + ig_b$sub2
## prec$sub2 <- rgamma(1, shape = shape_sub2, rate = rate_sub2)


## This is a Gibbs sampler v3 for longitudinal Bayesian semiparametric ridge.
## Update two block of parameters: variance and coefs

## Requirements: y, grp, Bmat (df), Xmat (df), Kmat, dim_block, ranef

## Algorithms paremeters: burn, size
## Extras: init (list of matrices with 'pop' and 'sub'), prior (see 'check_prior')


#' Bayesian ridge for longitudinal semiparemetric models
#'
#' This is a Gibbs sampler v3 for fitting Bayesian longitudinal semiparametric
#' models. The model is
#' y_{ij} = Bmat
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
#'   corresponding coefficients should be treated as random effects.
#' @param burn The number of samples to be burned before actual sampling.
#' @param size The number of samples to be obtained from the samples.
#' @param init unknown yet
#' @param prior unknown yet
#' @param prec unknown yet
#' 
#' @export
bayes_ridge_semi <- function(y, grp, Bmat, Xmat,
                             Kmat = NULL, dim_block = NULL, ranef = NULL,
                             burn = 0, size = 1000, init = NULL, prior = NULL,
                             prec = NULL, debug = TRUE) {
  
  grp <- purrr::map(grp, ~factor(.x, levels = unique(.x)))
  check_Kmat(Kmat)
  ## grp <- check_grp(grp)
  
  ## para: n_subs, subs_in_pop, pop_of_subs
  para <- get_para(grp, Bmat, Xmat, Kmat, dim_block, ranef)
  coef_samples <- get_empty_coef_samples(para, size)
  prec_samples <- get_empty_prec_samples(para, size)

  ## get a list where each element is a list of df split according to subjects
  datals <- split_data_subjects(Bmat, Xmat, y, grp$sub, debug)

  pls <- get_pls(y, Bmat$pop, Kmat)
  kcoef <- initialise_with_pls_v3(init, para, pls, datals, debug)

  ## some pre-calculation
  xX <- purrr::map(datals$Xmat, crossprod)
  prior_ls <- prior
  
  for (k in seq.int(-burn + 1, size)) {
    if (is.null(prec)) {
      kprec <- update_precs(kcoef, y, para, prior_ls)
    } else {
      kprec <- prec
    }
    kcoef <- update_coefs(datals, xX, para, kprec, debug)
  if (k > 0) {
      stopifnot(NROW(kcoef$theta) == para$dim_theta,
                NROW(kcoef$delta) == para$dim_delta,
                NROW(kcoef$beta) == para$dim_beta)
      coef_samples$population[, colnames(kcoef$theta), k] <- kcoef$theta
      coef_samples$subjects[, colnames(kcoef$delta), k] <- kcoef$delta
      coef_samples$beta[, k] <- kcoef$beta
      
      prec_samples$theta[, , k] <- kprec$theta
      prec_samples$delta[, , k] <- kprec$delta
      prec_samples$beta[, , k] <- kprec$beta
      prec_samples$eps[k] <- kprec$eps
    }
  }
  means <- purrr::map(coef_samples, ~rowMeans(.x, dims = length(dim(.x)) - 1))
  list(means = means, coef_samples = coef_samples, prec_samples = prec_samples)
}
