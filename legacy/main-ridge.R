## initialise theta with a penalised LS estimate, and delta with rnorms with small sd
initialise_with_pls <- function(init, para, grp, pls) {

  ## init: NULL or list(pop, sub); any of the element can be NULL
  n_terms_pop <- para$n_terms_pop
  n_terms_sub <- para$n_terms_sub
  n_pops <- para$n_pops
  n_subs <- para$n_subs

  if (is.null(init$pop)) {
    init$pop <- matrix(pls, n_terms_pop, n_pops,
                       dimnames = list(NULL, levels(grp$pop)))
  } else {
    if (all(dim(init$pop) == c(n_terms_pop, n_pops))) {
      init$pop <- as.matrix(init$pop)
      message("Population initial values supplied.")
      if (is.null(colnames(init$pop))) {
        colnames(init$pop) <- levels(grp$pop)
      } else {
        init$pop <- init$pop[, levels(grp$pop)]
      }
    } else {
      stop("Invalid dimension of population initial values.")
    }
  }
  if (is.null(init$sub)) {
    init$sub <- matrix(stats::rnorm(n_terms_sub * n_subs) * 0.01, n_terms_sub, n_subs,
                       dimnames = list(NULL, levels(grp$sub)))
  } else {
    if (all(dim(init$sub) == c(n_terms_sub, n_subs))) {
      init$sub <- as.matrix(init$sub)
      message("Subjects initial values supplied.")
      if (is.null(colnames(init$sub))) {
        colnames(init$sub) <- levels(grp$sub)
      } else {
        init$sub <- init$sub[, levels(grp$sub)]
      }
    } else {
      stop("Invalid dimension of subject initial values.")
    }
  }
  init
}

## calculates unnormalised log posterior for the subject specific model
## vector of coef_pop
## matrix of coef_sub
## list of prec: pop, sub1, sub2, eps
## list of contrib: pop, sub
## list of para: xKmat, rank_K, dim_sub1, y
## list of prior hyperparameters: ig_a, ig_b, iw_v, iw_lambda
logprior_sub <- function(coef_pop, coef_sub, prec, para, prior) {

  check_names <- prod(c("pop", "sub1", "sub2", "eps") %in% names(prec),
                      c("xKmat", "rank_K", "dim_sub1") %in% names(para),
                      c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                      c("pop", "sub2", "eps") %in% names(prior$ig_a),
                      c("pop", "sub2", "eps") %in% names(prior$ig_b))

  if (check_names == 0) {
    stop("Missing elements in the list.")
  }

  xKmat <- para$xKmat
  rank_K <- para$rank_K
  dim_sub1 <- para$dim_sub1
  ig_a <- prior$ig_a
  ig_b <- prior$ig_b
  iw_v <- prior$iw_v
  iw_lambda <- prior$iw_lambda          # assuming symmetric

  prec_sub <- block_diag(prec$sub1, diag(prec$sub2, NROW(coef_sub) - dim_sub1))

  ## inverse gamma
  lp_prec_pop <- (ig_a$pop + 1) * log(prec$pop) - ig_b$pop * prec$pop
  lp_prec_sub2 <- (ig_a$sub2 + 1) * log(prec$sub2) - ig_b$sub2 * prec$sub2
  lp_prec_eps <- (ig_a$eps + 1) * log(prec$eps) - ig_b$eps * prec$eps

  ## inverse Wishart
  lp_prec_sub1 <- 0.5 * (iw_v + dim_sub1 + 1) * determinant(prec$sub1)$modulus -
                                                                     0.5 * sum(diag(iw_lambda %*% prec$sub1))

  ## Gaussian
  lp_pop <- 0.5 * rank_K * log(prec$pop) -
    0.5 * prec$pop * crossprod(coef_pop, xKmat %*% coef_pop)
  lp_sub <- NCOL(coef_sub) / 2 * determinant(prec_sub)$modulus -
                                                     0.5 * sum(apply(coef_sub, 2, function(x) crossprod(x, prec_sub %*% x)))

  lp_prec_pop + lp_prec_sub1 + lp_prec_sub2 + lp_prec_eps + lp_pop + lp_sub
}

## calculates unnormalised log-likelihood for the subject specific model
## list of prec: eps
## resids vector
## n_samples
loglike_sub <- function(prec, resids, n_samples) {
  check_names <- prod(c("eps") %in% names(prec))

  if (check_names == 0) {
    stop("Missing elements in the list.")
  }

  n_samples / 2 * log(prec$eps) - 0.5 * prec$eps * crossprod(resids)
}

check_grp <- function(grp) {

  ## check grp field names
  if (!all(c("pop", "sub") %in% names(grp))) {
    stop("Missing grp variables.")
  }

  ## check grp lengths
  if (length(grp$pop) != length(grp$sub)) {
    stop("Unequal length of grp$pop and grp$sub")
  }

  ## sub factor must be nested within pop factor
  sublvl_in_pop <- tapply(grp$sub, grp$pop, unique, simplify = FALSE)

  if (length(unlist(sublvl_in_pop)) > length(unique(grp$sub))) {
    stop("Each subject can only appear in one population.")
  }
  grp <- lapply(grp, factor)
  lapply(grp, droplevels)
}

## Bmat is a list, Kmat is a matrix
check_Bmat <- function(Bmat, Kmat) {

  ## check grp field names
  if (!all(c("pop", "sub") %in% names(Bmat))) {
    stop("Missing Bmat variables.")
  }

  ## check Bmat lengths
  if (NROW(Bmat$pop) != NROW(Bmat$sub)) {
    stop("Unequal length of Bmat$pop and Bmat$sub.")
  }

  ## check col of Bmat
  if (NCOL(Bmat$pop) != NCOL(Kmat)) {
    stop("Inconsistant columns of Bmat and Kmat.")
  }
}

## prior is a list, dim_sub1 is an integer
check_prior <- function(prior, dim_sub1) {

  ## build a prior for covariance if not give
  if (is.null(prior)) {
    prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                  ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                  iw_v = dim_sub1 + 1,
                  iw_lambda = diag(dim_sub1)) # assuming symmetry
  }

  ## check prior field name
  check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                      c("pop", "sub2", "eps") %in% names(prior$ig_a),
                      c("pop", "sub2", "eps") %in% names(prior$ig_b))
  if (check_names == 0) {
    stop("Missing prior hyperparameters.")
  }
  prior
}

## Bmat, idx_sub and para are lists
calc_xB <- function(Bmat, idx_sub, para) {

  ## crossproduct of Bmats (pop x pop, sub x pop, sub x sub)
  n_subs <- para$n_subs
  n_terms_pop <- para$n_terms_pop
  n_terms_sub <- para$n_terms_sub

  xB_pop <- array(NA, c(n_terms_pop, n_terms_pop, n_subs),
                  list(NULL, NULL, names(idx_sub)))
  xB_subpop <- array(NA, c(n_terms_sub, n_terms_pop, n_subs),
                     list(NULL, NULL, names(idx_sub)))
  xB_sub <- array(NA, c(n_terms_sub, n_terms_sub, n_subs),
                  list(NULL, NULL, names(idx_sub)))

  for (i in names(idx_sub)) {
    xB_pop[, , i] <- crossprod(Bmat$pop[idx_sub[[i]], ])
    xB_subpop[, , i] <- crossprod(Bmat$sub[idx_sub[[i]], ], Bmat$pop[idx_sub[[i]], ])
    xB_sub[, , i] <- crossprod(Bmat$sub[idx_sub[[i]], ])
  }
  list(pop = xB_pop, subpop = xB_subpop, sub = xB_sub)
}

## y is a vector, Bmat, idx_sub and para are lists
calc_Bxy <- function(y, Bmat, idx_sub, para) {

  n_terms_pop <- para$n_terms_pop
  n_terms_sub <- para$n_terms_sub
  n_subs <- para$n_subs

  ## crossproduct of Bmat and y (Bmat_pop x y, Bmat_sub x y)
  Bxy_pop <- matrix(NA, n_terms_pop, n_subs, dimnames = list(NULL, names(idx_sub)))
  Bxy_sub <- matrix(NA, n_terms_sub, n_subs, dimnames = list(NULL, names(idx_sub)))

  for (i in names(idx_sub)) {
    Bxy_pop[, i] <- crossprod(Bmat$pop[idx_sub[[i]], ], y[idx_sub[[i]]])
    Bxy_sub[, i] <- crossprod(Bmat$sub[idx_sub[[i]], ], y[idx_sub[[i]]])
  }
  list(pop = Bxy_pop, sub = Bxy_sub)
}

## update precision
update_prec <- function(coef, resids, para, prior) {

  ## coef: list(pop = n_terms * 1 vec, sub = n_terms * n_subs mat)
  ## resid: n_samples * 1 vec
  ## para: list(Kmat, rank_K, dim_sub1, n_samples, n_subs, n_terms)
  ## prior: list(ig_a = list(pop, sub2, eps), ig_b = list(pop, sub2, eps),
  ##        iw_v, iw_lambda = dim_sub1 * dim_sub1 mat)

  Kmat <- para$Kmat
  rank_K <- para$rank_K
  n_samples <- para$n_samples
  n_pops <- para$n_pops
  n_subs <- para$n_subs
  n_terms_sub <- para$n_terms_sub
  dim_sub1 <- para$dim_sub1

  ig_a <- prior$ig_a                  # a list
  ig_b <- prior$ig_b                  # a list
  iw_v <- prior$iw_v                  # a number
  iw_lambda <- prior$iw_lambda        # a matrix, assuming symmetric

  prec <- list(pop = NA, eps = NA, sub1 = NA, sub2 = NA)

  ## update sigma^2_theta
  shape_pop <- 0.5 * n_pops * rank_K + ig_a$pop
  rate_pop <- 0.5 * sum((Kmat %*% coef$pop)^2) + ig_b$pop
  prec$pop <- stats::rgamma(1, shape = shape_pop, rate = rate_pop)

  ## update Sigma_dev1
  df_sub1 <- iw_v + n_subs
  scale_sub1 <- iw_lambda + tcrossprod(coef$sub[1:dim_sub1, , drop = FALSE])
  inv_scale_sub1 <- chol2inv(chol(scale_sub1))
  prec$sub1 <- stats::rWishart(1, df = df_sub1, Sigma = inv_scale_sub1)[, , 1]

  ## update sigma^2_dev2
  shape_sub2 <- 0.5 * n_subs * (n_terms_sub - dim_sub1) + ig_a$sub2
  rate_sub2 <- 0.5 * sum(coef$sub[-(1:dim_sub1), ]^2) + ig_b$sub2
  prec$sub2 <- stats::rgamma(1, shape = shape_sub2, rate = rate_sub2)
  
  ## update sigma^2_epsilon
  shape_eps <- 0.5 * n_samples + ig_a$eps
  rate_eps <- 0.5 * crossprod(resids) + ig_b$eps
  prec$eps <- stats::rgamma(1, shape = shape_eps, rate = rate_eps)

  prec
}

## update theta_l
update_theta_l <- function(xB_l, Bxy_l, xKmat, kprec, para) {

  n_subs_in_pop <- NCOL(Bxy_l$pop)
  n_terms_pop <- para$n_terms_pop

  ## xB_l (list of arrays of relevant xB)
  if (!all(c("pop", "sub", "subpop") %in% names(xB_l))) {
    stop("Missing xB_l variables.")
  }
  xB_sub_l <- xB_l$sub
  xB_pop_l <- xB_l$pop
  xB_subpop_l <- xB_l$subpop

  ## Bxy_l (list of matrices of relevant Bxy)
  if (!all(c("pop", "sub") %in% names(Bxy_l))) {
    stop("Missing xB_l variables.")
  }
  Bxy_pop_l <- Bxy_l$pop
  Bxy_sub_l <- Bxy_l$sub

  BMB <- array(NA, c(n_terms_pop, n_terms_pop, n_subs_in_pop))
  BMy <- matrix(NA, n_terms_pop, n_subs_in_pop)
  for (i in seq_len(n_subs_in_pop)) {
    Li <- xB_sub_l[, , i] + kprec$sub / kprec$eps
    inv_Li <- chol2inv(chol(Li))
    BBLBB_i <- crossprod(xB_subpop_l[, , i], inv_Li %*% xB_subpop_l[, , i])
    BMB[, , i] <- kprec$eps * (xB_pop_l[, , i] - BBLBB_i)
    BBLBy_i <- crossprod(xB_subpop_l[, , i], inv_Li %*% Bxy_sub_l[, i])
    BMy[, i] <- kprec$eps * (Bxy_pop_l[, i] - BBLBy_i)
  }
  Phi <- kprec$pop * xKmat + rowSums(BMB, dims = 2)
  inv_Phi <- chol2inv(chol(Phi))
  t(mvtnorm::rmvnorm(1, inv_Phi %*% rowSums(BMy), inv_Phi))
}


## update delta_i
update_delta_i <- function(Bmat_sub_i, xB_sub_i, y_star_i, kprec, para) {
  M_sub <- chol2inv(chol(xB_sub_i + kprec$sub / kprec$eps))
  mu <- M_sub %*% crossprod(Bmat_sub_i, y_star_i)
  sig <- M_sub / kprec$eps
  t(mvtnorm::rmvnorm(1, mu, sig))
}

## This is a Gibbs sampler v2 for longitudinal Bayesian ridge; see thesis.
## Update two block of parameters: variance and coefs
## Requirements: Bmat (list), y, grp (list), Kmat, dim_sub1
## Algorithms paremeters: burn, size
## Extras: init (list of matrices with 'pop' and 'sub'), prior (see 'check_prior')
#' @export
bayes_ridge_sub_v2 <- function(y, grp, Bmat, Kmat, dim_sub1, burn, size,
                               init = NULL, prior = NULL) {

  grp <- check_grp(grp)
  prior <- check_prior(prior, dim_sub1)
  check_Bmat(Bmat, Kmat)

  ## some precalculation
  rank_K <- NROW(Kmat)
  n_terms_pop <- NCOL(Bmat$pop)
  n_terms_sub <- NCOL(Bmat$sub)
  n_samples <- NROW(Bmat$sub)
  n_pops <- length(unique(grp$pop))
  n_subs <- length(unique(grp$sub))

  para <- list(Kmat = Kmat, dim_sub1 = dim_sub1, rank_K = rank_K,
               n_terms_pop = n_terms_pop, n_terms_sub = n_terms_sub,
               n_samples = n_samples, n_pops = n_pops, n_subs = n_subs)

  subs_in_pop <- tapply(as.character(grp$sub), grp$pop, unique, simplify = FALSE)
  idx_pop <- tapply(seq_len(n_samples), grp$pop, function(x) x, simplify = FALSE)
  idx_sub <- tapply(seq_len(n_samples), grp$sub, function(x) x, simplify = FALSE)

  ## some crossproducts precalculation
  xB <- calc_xB(Bmat, idx_sub, para)      # pop x pop, sub x pop, sub x sub
  Bxy <- calc_Bxy(y, Bmat, idx_sub, para) # Bmat_pop x y, Bmat_sub x y
  xKmat <- crossprod(Kmat)

  ## initialise theta with a penalised LS estimate, and delta with rnorms with
  ## small sd
  init <- initialise_with_pls(init, para, grp, get_pls(y, Bmat$pop, Kmat))
  kcoef <- list(pop = init$pop, sub = init$sub)

  ## initialise prediction contribution by population coefs and subjects deviations
  kcontrib_pop <- rep(NA, n_samples)
  for (l in levels(grp$pop)) {
    kcontrib_pop[idx_pop[[l]]] <- Bmat$pop[idx_pop[[l]], ] %*% kcoef$pop[, l]
  }
  kcontrib_sub <- rep(NA, n_samples)
  for (i in levels(grp$sub)) {
    kcontrib_sub[idx_sub[[i]]] <- Bmat$sub[idx_sub[[i]], ] %*% kcoef$sub[, i]
  }

  ## initialise the output list
  samples <- initialise_samples(para, size, dim_sub1, grp)

  for (k in seq.int(-burn + 1, size)) {
    ## update precisions
    kresids <- y - kcontrib_pop - kcontrib_sub
    kprec <- update_prec(kcoef, kresids, para, prior)
    kprec$sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms_sub - dim_sub1))

    ## update theta
    for (l in levels(grp$pop)) {
      xB_l <- list(pop = xB$pop[, , subs_in_pop[[l]]],
                   sub = xB$sub[, , subs_in_pop[[l]]],
                   subpop = xB$subpop[, , subs_in_pop[[l]]])
      Bxy_l <- list(pop = Bxy$pop[, subs_in_pop[[l]]],
                    sub = Bxy$sub[, subs_in_pop[[l]]])
      kcoef$pop[, l] <- update_theta_l(xB_l, Bxy_l, xKmat, kprec, para)
      kcontrib_pop[idx_pop[[l]]] <- Bmat$pop[idx_pop[[l]], ] %*% kcoef$pop[, l]
    }

    ## update delta
    for (i in levels(grp$sub)) {
      Bmat_sub_i <- Bmat$sub[idx_sub[[i]], ]
      xB_sub_i <- xB$sub[, , i]
      y_star_i <- y[idx_sub[[i]]] - kcontrib_pop[idx_sub[[i]]]
      kcoef$sub[, i] <- update_delta_i(Bmat_sub_i, xB_sub_i, y_star_i, kprec, para)
      kcontrib_sub[idx_sub[[i]]] <- Bmat$sub[idx_sub[[i]], ] %*% kcoef$sub[, i]
    }

    ## print progress
    if (k %% 1000 == 0) {
      message(k, " samples generated.")
    }

    ## fix from here on
    ## store samples after burn-in iterations
    if (k > 0) {
      samples$precision$pop[k] <- kprec$pop
      samples$precision$sub1[, , k] <- kprec$sub1
      samples$precision$sub2[k] <- kprec$sub2
      samples$precision$eps[k] <- kprec$eps
      samples$population[, , k] <- kcoef$pop
      samples$subjects[, , k] <- kcoef$sub

      ## ## calculate unnormalised log-likelihood
      ## kresids <- y - kcontrib_pop - kcontrib_sub
      ## samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

      ## ## calculate unnormalised log-posterior
      ## para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
      ## samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
      ##     samples$ll[k]
    }
  }
  means <- list(population = rowMeans(samples$population, dims = 2),
                subjects = rowMeans(samples$subjects, dims = 2))
  list(means = means, samples = samples)
}

## This is a Gibbs sampler for constrained longitudinal Bayesian ridge, with an
## HMC step to generate pop and sub coefs; see thesis.
## A * theta >= 0 ; A * (theta + delta_i) >= 0
## The vector of constriant has to be zero.
## Assumption: Full row rank for Amat.
## Requirements: Bmat, y, grp, Kmat, dim_sub1, Amat
## Algorithms paremeters: burn, size

#' @export
bayes_ridge_cons_sub_v2 <- function(y, grp, Bmat, Kmat, dim_sub1, Amat, burn, size,
                                    init = NULL, prior = NULL, prec = NULL) {

  ## hyperparemeters for priors
  if (is.null(prior)) {
    prior <- list(ig_a = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                  ig_b = list(pop = 0.001, sub2 = 0.001, eps = 0.001),
                  iw_v = dim_sub1 + 1,
                  iw_lambda = diag(dim_sub1)) # assuming symmetry
  }

  check_names <- prod(c("ig_a", "ig_b", "iw_v", "iw_lambda") %in% names(prior),
                      c("pop", "sub2", "eps") %in% names(prior$ig_a),
                      c("pop", "sub2", "eps") %in% names(prior$ig_b))
  if (check_names == 0) {
    stop("Missing prior hyperparameters.")
  }

  if (NCOL(Amat) != NCOL(Bmat)) {
    stop("Dims of Amat and Bmat inconsistent.")
  }

  ## clean up grp variable
  if (is.factor(grp)) {
    grp <- droplevels(grp)
  } else {
    grp <- factor(grp, levels = unique(grp))
  }

  ## some precalculation
  rank_K <- NROW(Kmat)
  n_terms <- NCOL(Bmat)
  n_samples <- NROW(Bmat)
  n_subs <- length(unique(grp))
  idx <- tapply(seq_len(n_samples), grp, function(x) x, simplify = FALSE)

  ## construct full constraint matrix
  fAmat_left <- t(matrix(t(Amat), NCOL(Amat), NROW(Amat) * (n_subs + 1)))
  fAmat_topright <- matrix(0, NROW(Amat), n_terms * n_subs)
  fAmat_btmright <- block_diag(Amat, size = n_subs)
  fAmat <- cbind(fAmat_left, rbind(fAmat_topright, fAmat_btmright))
  lower <- rep(0, NROW(Amat))
  flower <- rep(lower, n_subs + 1)

  ## construct full design matrix
  fBmat_right <- do.call(block_diag, lapply(idx, function(x) Bmat[x, ]))
  fBmat <- cbind(Bmat, fBmat_right)

  ## some crossproducts precalculation
  xfBmat <- crossprod(fBmat)
  fBxy <- crossprod(fBmat, y)
  xKmat <- crossprod(Kmat)

  ## initialise the output list
  samples <- list(population = matrix(NA, n_terms, size),
                  subjects = array(NA, c(n_terms, n_subs, size),
                                   dimnames = list(NULL, levels(grp))),
                  precision = list(pop = rep(NA, size),
                                   sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                   sub2 = rep(NA, size),
                                   eps = rep(NA, size)),
                  lp = rep(NA, size),
                  ll = rep(NA, size))

  ## initialise theta and delta with tnorms with small sd
  init <- initialise_with_Amat(init, n_terms, grp, Amat)
  kcoef <- list(pop = init$pop, sub = init$sub)
  fkcoef <- c(kcoef$pop, kcoef$sub)

  ## initialise prediction contribution by population coefs and subjects deviations
  kcontrib_pop <- Bmat %*% kcoef$pop
  kcontrib_sub <- rep(NA, n_samples)
  for (i in levels(grp)) {
    kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
  }

  for (k in seq.int(-burn + 1, size)) {
    ## update precisions
    para <- list(Kmat = Kmat, rank_K = rank_K, dim_sub1 = dim_sub1,
                 n_samples = n_samples, n_subs = n_subs, n_terms = n_terms)
    kresids <- y - kcontrib_pop - kcontrib_sub
    kprec <- update_prec(kcoef, kresids, para, prior)

    ## fix precision if provided
    if (!is.null(prec)) {
      kprec$pop <- prec$pop
      kprec$sub1 <- prec$sub1
      kprec$sub2 <- prec$sub2
    }

    ## construct Sigma
    kprec_sub <- block_diag(kprec$sub1, diag(kprec$sub2, n_terms - dim_sub1))
    fkprec_all <- block_diag(kprec_sub, size = n_subs + 1)
    fkprec_all[1:n_terms, 1:n_terms] <- kprec$pop * xKmat

    ## update theta and delta
    M <- chol2inv(chol(xfBmat + fkprec_all / kprec$eps))
    mu <- M %*% fBxy
    sig <- 1 / kprec$eps * M
    fkcoef <- t(tnorm::rmvtnorm(1, mu, sig, fkcoef, fAmat, -flower, 10))
    ## fkcoef <- t(mvtnorm::rmvnorm(1, mu, sig))  # no problem in unconstrained case
    kcoef$pop <- fkcoef[1:n_terms]
    kcoef$sub <- matrix(fkcoef[-(1:n_terms)], n_terms, n_subs,
                        dimnames = list(NULL, levels(grp)))
    kcontrib_pop <- Bmat %*% kcoef$pop
    for (i in levels(grp)) {
      kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef$sub[, i]
    }

    ## print progress
    if (k %% 1000 == 0) {
      message(k, " samples generated.")
    }

    ## store samples after burn-in iterations
    if (k > 0) {
      samples$precision$pop[k] <- kprec$pop
      samples$precision$sub1[, , k] <- kprec$sub1
      samples$precision$sub2[k] <- kprec$sub2
      samples$precision$eps[k] <- kprec$eps
      samples$population[, k] <- kcoef$pop
      samples$subjects[, , k] <- kcoef$sub

      ## calculate unnormalised log-likelihood
      kresids <- y - kcontrib_pop - kcontrib_sub
      samples$ll[k] <- loglike_sub(kprec, kresids, n_samples)

      ## calculate unnormalised log-posterior
      para <- list(xKmat = xKmat, rank_K = rank_K, dim_sub1 = dim_sub1)
      samples$lp[k] <- logprior_sub(kcoef$pop, kcoef$sub, kprec, para, prior) +
        samples$ll[k]
    }
  }

  means <- list(population = rowMeans(samples$population),
                subjects = rowMeans(samples$subjects, dims = 2))
  list(means = means, samples = samples)
}

