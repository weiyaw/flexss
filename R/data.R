get_simdata2 <- function(seed = 1) {
  ## a quadratic spline with 3 interior knots
  set.seed(seed)
  x <- 1:20 + rnorm(20, sd = 0.1) # 20 samples for each subject
  K <- 3
  deg <- 2
  n_pops <- 3
  n_subs <- 10
  ## 1st group = 3 subs, 2nd group = 3 subs, 3rd group = 4 subs, 
  n_subs_in_pop <- list('1' = 3, '2' = 3, '3' = 4)
  grp <- list(pop = unlist(purrr::imap(n_subs_in_pop, ~rep(.y, .x * length(x))), F, F),
              sub = as.character(rep(1:10, each = length(x))))

  Bmat <- get_design_bs(x, K, deg)$design
  
  dim_theta <- K + deg + 1
  dim_delta <- K + deg + 1

  ## theta: 1:dim_theta + std Gaussian
  ## delta: std Gaussian
  theta <- matrix(1:dim_theta + rnorm(n_pops * dim_theta), dim_theta, n_pops,
                  dimnames = list(NULL, seq(1, n_pops)))
  delta <- matrix(rnorm(n_subs * dim_delta) * 0.5, dim_delta, n_subs)
  
  f <- unlist(purrr::imap(n_subs_in_pop, ~rep(Bmat %*% theta[, .y], .x)), F, F)
  g <- c(Bmat %*% delta)
  y <- f + g + rnorm(length(f), sd = 0.2)

  data.frame(y = y, x = rep(x, n_subs),  sub = grp$sub, pop = grp$pop, truth = f + g)
}


get_simdata3 <- function(seed = 1) {
  ## a quadratic spline with 3 interior knots
  set.seed(seed)
  x <- 1:20 + rnorm(20, sd = 0.1) # 20 samples for each subject
  K <- 3
  deg <- 2
  n_pops <- 3
  n_subs <- 10
  ## 1st group = 3 subs, 2nd group = 3 subs, 3rd group = 4 subs, 
  n_subs_in_pop <- list('1' = 3, '2' = 3, '3' = 4)
  grp <- list(pop = unlist(purrr::imap(n_subs_in_pop, ~rep(.y, .x * length(x))), F, F),
              sub = as.character(rep(1:10, each = length(x))))

  Bmat <- get_design_bs(x, K, deg)$design

  ## effect1: sub 1, 3, 5, 7, 9 has level 1, the rest level 0
  effect1 <- as.factor(as.numeric(grp$sub) %% 2)
  ## effect2: sub 1, 4, 7, 10 has level 1, sub 2, 5, 8 has level 2, the rest level 0
  effect2 <- as.factor(as.numeric(grp$sub) %% 3)
  dim_theta <- K + deg + 1
  dim_delta <- K + deg + 1

  ## theta: 1:dim_theta + std Gaussian
  ## delta: std Gaussian
  theta <- matrix(1:dim_theta + rnorm(n_pops * dim_theta), dim_theta, n_pops,
                  dimnames = list(NULL, seq(1, n_pops)))
  delta <- matrix(rnorm(n_subs * dim_delta) * 0.5, dim_delta, n_subs)
  
  f <- unlist(purrr::imap(n_subs_in_pop, ~rep(Bmat %*% theta[, .y], .x)), F, F)
  g <- c(Bmat %*% delta)
  y <- f + g + rnorm(length(f), sd = 0.2)

  data.frame(y = y, x = rep(x, n_subs),  sub = grp$sub, pop = grp$pop, truth = f + g,
             effect1 = effect1, effect2 = effect2)
}

