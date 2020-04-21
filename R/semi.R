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
    Nmat_i <- crossprod(Bmat_sub_i, PhixBs_i) + kprec$sub / kprec$eps
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
      purrr::reduce(`+`, .init = kprec$pop / kprec$eps) # this step gives Q_l
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
calc_f <- function(Bmat_pop, theta, pop_of_subs) {
  purrr::map2(Bmat_pop, pop_of_subs, ~.x %*% theta[[.y]])
}

## return: list of length n_sub
calc_g <- function(Bmat_sub, delta) {
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

## return a 1-row matrix
update_beta <- function(Lmat, Mmat, y, f, g, kprec) {
  Lmat_sum <- purrr::reduce(Lmat, `+`)
  Lmat_sum_inv <- chol2inv(chol(Lmat_sum))
  Mmat_sum <- purrr::reduce(Mmat, `+`)
  mu_beta <- tcrossprod(Lmat_sum_inv, Mmat_sum)
  Sig_beta <- Lmat_sum_inv / kprec$eps
  mvtnorm::rmvnorm(1, mu_beta, Sig_beta)
  ## notability notes typo: there should be an epsilon in the cov matrix
}

## para: n_subs, subs_in_pop, pop_of_subs
update_coefs <- function(Bmat_pop, Bmat_sub, Xmat, y, xX, para, kprec, debug = F) {
  if (debug) {
    stopifnot(names(Bmat_pop) == names(Bmat_sub),
              names(Bmat_pop) == names(Xmat),
              names(Bmat_pop) == names(y),
              names(Bmat_pop) == names(xX),
              names(Bmat_pop) == names(para$pop_of_subs))
  }

  ## some pre-calculation
  Lmat <- calc_Lmat(xX, kprec, para$n_subs)
  Phi <- calc_Phi(Xmat, Lmat)
  PhixBs <- calc_PhixBs(Phi, Bmat_sub)
  Nmat_inv <- calc_Nmat_inv(Bmat_sub, PhixBs, kprec)

  ## update theta
  Psy <- calc_Psy(Phi, PhixBs, Nmat_inv)
  PsyxB <- calc_PsyxB(Psy, Bmat_pop)
  theta <- update_theta(Bmat_pop, PsyxB, y, para$subs_in_pop, kprec)

  ## update delta
  f <- calc_f(Bmat_pop, theta, para$pop_of_subs)
  delta <- update_delta(Nmat_inv, PhixBs, y, f, kprec)
  
  ## update beta
  g <- calc_g(Bmat_sub, delta)
  Mmat <- calc_Mmat(y, f, g, Xmat)
  beta <- update_beta(Lmat, Mmat, y, f, g, kprec)
  
  ## theta & delta are lists of row vectors
  list(theta = do.call(cbind, theta) %>% magrittr::set_colnames(names(theta)),
       delta = do.call(cbind, delta) %>% magrittr::set_colnames(names(delta)),
       beta = beta)
}
  
get_empty_coef_samples <- function(para, size) {
  pop_array <- array(NA, c(para$n_terms_pop, para$n_pops, size),
                     dimnames = list(NULL, names(para$subs_in_pop)))
  sub_array <- array(NA, c(para$n_terms_sub, para$n_subs, size),
                     dimnames = list(NULL, names(para$pop_of_subs)))
  beta_matrix <- matrix(NA, para$n_terms_beta, size)
  list(population = pop_array,
       subjects = sub_array,
       beta = beta_matrix,
       lp = rep(NA, size),
       ll = rep(NA, size))
}

get_para <- function(grp, Bmat, Xmat) {
  para <- list(n_pops = length(unique(grp$pop)),
               n_subs = length(unique(grp$sub)),
               n_terms_pop = NCOL(Bmat$pop),
               n_terms_sub = NCOL(Bmat$sub),
               n_terms_beta = NCOL(Xmat))
  para$subs_in_pop <- tapply(grp$sub, grp$pop,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  para$pop_of_subs <- tapply(grp$pop, grp$sub,
                             function(x) as.character(unique(x)),
                             simplify = FALSE)
  para
}

## This is a Gibbs sampler v3 for longitudinal Bayesian semiparametric ridge.
## Update two block of parameters: variance and coefs
## Requirements: Bmat (list), y, grp (list), Kmat, dim_sub1
## Algorithms paremeters: burn, size
## Extras: init (list of matrices with 'pop' and 'sub'), prior (see 'check_prior')
bayes_ridge_semi <- function(y, grp, Bmat, Xmat, Kmat = NULL, dim_sub1 = NULL,
                             burn = 0, size = 1000, init = NULL, prior = NULL,
                             prec = NULL) {
  
  grp <- purrr::map(grp, ~factor(.x, levels = unique(.x)))
  ## grp <- check_grp(grp)
  
  ## para: n_subs, subs_in_pop, pop_of_subs
  para <- get_para(grp, Bmat, Xmat)

  Bmat_pop <- split.data.frame(Bmat$pop, grp$sub)
  Bmat_sub <- split.data.frame(Bmat$sub, grp$sub)
  Xmat <- split.data.frame(Xmat, grp$sub)
  y <- split.default(y, grp$sub)
  xX <- purrr::map(Xmat, crossprod)
  
  coef_samples <- get_empty_coef_samples(para, size)

  kprec <- prec
  if (is.null(kprec)) {
    ## kprec <- list()
  }
  
  for (k in seq.int(-burn + 1, size)) {
    kcoef <- update_coefs(Bmat_pop, Bmat_sub, Xmat, y, xX, para, kprec, debug = TRUE)
    if (k > 0) {
      coef_samples$population[, colnames(kcoef$theta), k] <- kcoef$theta
      coef_samples$subjects[, colnames(kcoef$delta), k] <- kcoef$delta
      coef_samples$beta[, k] <- t(kcoef$beta)
    }
  }
  means <- list(population = rowMeans(coef_samples$population, dims = 2),
                subjects = rowMeans(coef_samples$subjects, dims = 2))
  list(means = means, samples = coef_samples)
}



