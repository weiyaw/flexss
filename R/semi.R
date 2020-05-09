
## the rows of Kmat must be independent
check_Kmat <- function(Kmat) {
  if (abs(det(tcrossprod(Kmat))) < 1e-10) {
    warning('The rows of Kmat may not be independent.')
  }
}

check_prec_beta <- function(prec_beta, dim_beta) {
  if (!is.null(prec$beta)) {
    message('prec$beta specified.')
    stopifnot(is.matrix(prec$beta),
              dim(prec$beta) == para$dim_beta,
              is_precision(prec$beta))
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
  update_yhat(f_init, g_init, datals$Xmat, beta)
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

calc_xX <- function(Xmat) {
  if (is.null(Xmat)) {
    NULL
  } else {
    purrr::map(Xmat, crossprod)
  }
}

## return: list of length n_sub
calc_Lmat <- function(xX, kprec, n_subs) {
  if (is.null(xX)) {
    NULL
  } else if (all(kprec$beta == 0)) {
    xX
  } else {
    purrr::map(xX, ~.x + kprec$beta / (n_subs * kprec$eps))
  }
}

## return: list of length n_sub
calc_Mmat <- function(y, f, g, Xmat) {
  if (is.null(Xmat)) {
    NULL
  } else {
    purrr::pmap(list(y, f, g, Xmat), ~crossprod(..1 - ..2 - ..3, ..4))
  }
}

## return: list of length n_sub
## if Xmat is NULL (i.e. no beta term), y must be provided
calc_Phi <- function(Xmat, Lmat, y = NULL) {
  if (is.null(Xmat)) {
    stopifnot(!is.null(y))
    purrr::map(y, ~diag(length(.x)))
  } else {
    calc_Phi_i <- function(Xmat_i, Lmat_i) {
      dim_Phi <- NROW(Xmat_i)
      Lmat_i_inv <- chol2inv(chol(Lmat_i)) # L_i > 0
      diag(dim_Phi) - tcrossprod(Xmat_i %*% Lmat_i_inv, Xmat_i)
    }
    purrr::map2(Xmat, Lmat, calc_Phi_i)
  }
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
  if (is.null(Lmat) || is.null(Mmat)) {
    NULL
  } else {
    Lmat_sum <- purrr::reduce(Lmat, `+`)
    Lmat_sum_inv <- chol2inv(chol(Lmat_sum))
    Mmat_sum <- purrr::reduce(Mmat, `+`)
    mu_beta <- tcrossprod(Lmat_sum_inv, Mmat_sum)
    Sig_beta <- Lmat_sum_inv / kprec$eps
    t(mvtnorm::rmvnorm(1, mu_beta, Sig_beta))
  }
}
## beta: a vector / column vector of length = NCOL(Xmat)
## the rest of args: list of column vectors/matrix, length = n_subs
## return a column vector (not numeric vector)
update_yhat <- function(f, g, Xmat, beta) {
  ## stopifnot(names(f) == names(g))
  ## if (is.null(Xmat)) {
  ##   unlist(f) + unlist(g)
  ## } else {
  ##   stopifnot(names(f) == names(Xmat))
  ##   tall_Xmat <- do.call(rbind, Xmat)
  ##   unlist(f) + unlist(g) + tall_Xmat %*% beta
  ## }

  ## safer to produce a list, or risk messing up the subject names
  stopifnot(names(f) == names(g))
  if (is.null(Xmat)) {
    purrr::map2(f, g, `+`)
  } else {
    stopifnot(names(f) == names(Xmat))
    purrr::pmap(list(f, g, Xmat), ~..1 + ..2 + ..3 %*% beta)
  }
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
  Phi <- calc_Phi(datals$Xmat, Lmat, datals$y)
  PhixBs <- purrr::map2(Phi, datals$Bmat_sub, `%*%`) # Phi_i is symmetrical
  Nmat_inv <- calc_Nmat_inv(datals$Bmat_sub, PhixBs, kprec)

  ## update theta
  Psy <- calc_Psy(Phi, PhixBs, Nmat_inv)
  PsyxB <- purrr::map2(Psy, datals$Bmat_pop, `%*%`) # Psy_i is symmetrical
  theta <- update_theta(datals$Bmat_pop, PsyxB, datals$y, para$subs_in_pop, kprec)

  ## update delta
  f <- purrr::map2(datals$Bmat_pop, para$pop_of_subs, ~.x %*% theta[[.y]])
  delta <- update_delta(Nmat_inv, PhixBs, datals$y, f, kprec)
  
  ## update beta (a column vector)
  g <- purrr::map2(datals$Bmat_sub, delta, `%*%`)
  Mmat <- calc_Mmat(datals$y, f, g, datals$Xmat)
  beta <- update_beta(Lmat, Mmat, datals$y, f, g, kprec)
  
  ## calc prediction
  yhat <- update_yhat(f, g, datals$Xmat, beta)

  ## theta & delta are matrices, each column represents one population/subject
  ## beta: a column vector
  ## yhat: a list of column vectors, length = n_subs
  list(theta = do.call(cbind, theta) %>% `colnames<-`(names(theta)),
       delta = do.call(cbind, delta) %>% `colnames<-`(names(delta)),
       beta = beta,
       yhat = yhat)
}

get_coef_container <- function(para, size) {
  theta_array <- array(NA, c(para$dim_theta, para$n_pops, size),
                       dimnames = list(NULL, names(para$subs_in_pop)))
  delta_array <- array(NA, c(para$dim_delta, para$n_subs, size),
                       dimnames = list(NULL, names(para$pop_of_subs)))
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
  para$subs_in_pop <- tapply(grp$sub, grp$pop,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  para$pop_of_subs <- tapply(grp$pop, grp$sub,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  
  if (is.null(Xmat)) {
    ## for ranef to be NULL if no random/fixed effects
    para[c('dim_beta', 'ranef', 'beta_names')] <- NULL
  } else {
    if (!is.null(prec_beta)) {
      message("ranef overriden by prec$beta") 
      para$ranef <- which(diag(prec$beta) > 0)
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
  resids <- unlist(y) - unlist(yhat)
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


## split data into chucks according to the subject index
## arguments are same as bayes_ridge_semi()
split_data_subjects <- function(Bmat, Xmat, y, grp_sub, debug) {
  datals <- list(Bmat_pop = split.data.frame(Bmat$pop, grp_sub),
                 Bmat_sub = split.data.frame(Bmat$sub, grp_sub),
                 y = split.default(y, grp_sub))

  if (is.null(Xmat)) {
    datals$Xmat <- NULL
  } else {
    datals$Xmat <- split.data.frame(Xmat, grp_sub)
  }

  if (debug) {
    stopifnot(names(datals$Bmat_pop) == names(datals$Bmat_sub),
              names(datals$Bmat_pop) == names(datals$Xmat),
              names(datals$Bmat_pop) == names(datals$y))
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
  
  grp <- purrr::map(grp, as.factor)
  check_Kmat(Kmat)
  para <- get_para(grp, Bmat, Xmat, Kmat, dim_block, ranef, prec$beta)
  check_prec(prec, para)
  ## grp <- check_grp(grp)
  
  ## para: n_subs, subs_in_pop, pop_of_subs
  coef_samples <- get_coef_container(para, size)
  prec_samples <- get_prec_container(para, size, prec)

  ## get a list where each element is a list of df split according to subjects
  datals <- split_data_subjects(Bmat, Xmat, y, grp$sub, debug)

  pls <- get_pls(y, Bmat$pop, Kmat)
  kcoef <- initialise_with_pls_v3(init, para, pls, datals, debug)

  ## some pre-calculation
  xX <- calc_xX(datals$Xmat)
  prior_ls <- prior
  
  for (k in seq.int(-burn + 1, size)) {
    kprec <- update_precs(kcoef, datals$y, para, prior_ls, prec)
    kcoef <- update_coefs(datals, xX, para, kprec, debug)

    ## print progress
    if (k %% 1000 == 0) {
      message(k, " samples generated.")
    }

    if (k > 0) {
      coef_samples$theta[, colnames(kcoef$theta), k] <- kcoef$theta
      coef_samples$delta[, colnames(kcoef$delta), k] <- kcoef$delta
      coef_samples$beta[, k] <- kcoef$beta
      
      prec_samples$theta[, , k] <- kprec$theta
      prec_samples$delta[, , k] <- kprec$delta
      prec_samples$beta[, , k] <- kprec$beta
      prec_samples$eps[k] <- kprec$eps
    }
  }
  
  list(samples = list(coef = coef_samples,
                      prec = prec_samples),
       means = list(coef = purrr::map(coef_samples, pmean),
                    prec = purrr::map(prec_samples, pmean)))
}
