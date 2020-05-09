################################################################
##########################            ##########################
##########################    PLOT    ##########################
##########################            ##########################
################################################################

## extract coefs from the model
extract_coef <- function(model, plot_type) {

  if (plot_type == "mean") {
    ## extract the posterior means
    coef_pop <- model$means$population
    dev_sub <- model$means$subjects
  } else if (plot_type == "map") {
    ## extract the map
    model_map <- get_max(model, "map")
    coef_pop <- model_map$population
    dev_sub <- model_map$subjects
  } else if (plot_type == "mle") {
    ## extract the mle (maximum likelihood)
    model_mle <- get_max(model, "mle")
    coef_pop <- model_mle$population
    dev_sub <- model_mle$subjects
  } else {
    stop("Unknown plotting type.")
  }

  list(pop = coef_pop, sub = dev_sub)
}

## calculate model matrix given x and a basis information
get_model_mat <- function(x, binfo) {

  knots <- binfo$knots
  names(knots) <- NULL
  deg <- binfo$degree

  if (binfo$type == "tpf") {
    model_mat <- get_design_tpf(x, knots, deg)$design
  } else if (binfo$type == "bs-ridge") {
    model_mat <- splines::splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
  } else if (binfo$type == "bs") {
    model_mat <- splines::splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
    model_mat <- model_mat %*% binfo$trans_mat
  } else {
    stop("Unknown type of basis.")
  }
  model_mat
}

get_plot_which <- function(plot_which, lvl_pop) {

  ## Plot all population if 'plot_which' is missing
  if (is.null(plot_which)) {
    plot_which <- lvl_pop
  } else if (!is.vector(plot_which)) {
    stop("'plot_which' must be a vector.")
  } else if (!all(plot_which %in% lvl_pop)) {
    stop("Invalid population levels in 'plot_which'.")
  }
  plot_which
}

get_limits <- function(limits, x) {

  ## Use the range of the predictor if 'limits' is missing
  if (is.null(limits)) {
    limits <- range(x)
  } else if (!(is.vector(limits) && length(limits) == 2)) {
    stop("'limits' must be a vector of length 2.")
  }
  limits
}

get_plot_x <- function(model, limits) {
  fine <- 200                         # how fine the plot_x should be?
  plot_x <- c(seq(min(limits), max(limits), length.out = fine),
              model$pop$basis$knots, model$sub$basis$knots)
  sort(unique(plot_x))
}

## get a dataset for ribbon plotting
get_plotdat_thin <- function(samples, model_mat, plot_x, thin = 10) {
  ## thin the sample by 10 by default
  thin_idx <- seq(thin, NCOL(samples), by = thin)
  thin_y <- model_mat %*% samples[, thin_idx]
  colnames(thin_y) <- thin_idx
  prob <- c(0.25, 0.75)
  prob_outer <- c(0.05, 0.95)
  apply(thin_y, 1, function(x) quantile(x, c(prob, prob_outer))) %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(x = plot_x)
}



## Plot the population and subject curves (suitable for models with multiple populations)
## 'plot_which' takes a vector of names (for models with multiple populations)
## type can be either of "mean", "mle" or "map"
## If shade is true, thin the population samples by 10 and plot the thinned samples.

## Requires: model$basis$knots (including extrema), model$basis$type,
##   model$basis$degree, model$info$lvl_pop, model$data, model$mean,
##   model$samples$population (if shade is true)

#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes_ geom_point geom_line geom_ribbon
plot_spline <- function(model, limits = NULL, plot_which = NULL, plot_type = "mean",
                        shade = FALSE, silent = FALSE) {

  data <- model$data
  subs_in_pop <- tapply(as.character(data$sub), data$pop, unique, simplify = FALSE)

  limits <- get_limits(limits, data$x)                    # limit of x
  plot_x <- get_plot_x(model, limits)                     # range of x to plot
  model_mat_pop <- get_model_mat(plot_x, model$basis$pop) # model matrix for the pop curves
  model_mat_sub <- get_model_mat(plot_x, model$basis$sub) # model matrix for the sub curves

  plot_which <- get_plot_which(plot_which, levels(data$pop)) # plot which populations?
  coefs <- extract_coef(model, plot_type)           # extract coefs and devs


  ggls <- list()
  for (l in plot_which) {
    ## shading
    if (shade) {
      plotdat_thin <- get_plotdat_thin(model$samples$population[, l, ],
                                       model_mat_pop, plot_x)
      ggls[[l]]$rib90 <- geom_ribbon(aes_(~x, ymin = ~`5%`, ymax = ~`95%`),
                                     plotdat_thin, fill = "grey85", alpha = 0.7)
      ggls[[l]]$rib50 <- geom_ribbon(aes_(~x, ymin = ~`25%`, ymax = ~`75%`),
                                     plotdat_thin, fill = "grey65", alpha = 0.5)
    }

    ## data
    ggls[[l]]$data <- geom_point(aes_(~x, ~y, col = ~sub), data[data$pop == l, ])

    ## curves
    subs_in_l <- subs_in_pop[[l]]
    plot_y_pop <- as.numeric(model_mat_pop %*% coefs$pop[, l])
    plot_y_sub <- plot_y_pop + (model_mat_sub %*% coefs$sub[, subs_in_l])
    plotdat_pop <- tibble::tibble(x = plot_x, y = plot_y_pop)
    plotdat_sub <- tibble::as_tibble(plot_y_sub) %>%
      dplyr::mutate(x = plot_x) %>%
      tidyr::gather("sub", "y", -"x")
    ggls[[l]]$sub <- geom_line(aes_(~x, ~y, col = ~sub, group = ~sub), plotdat_sub)
    ggls[[l]]$pop <- geom_line(aes_(~x, ~y), data = plotdat_pop)

    if (!silent) {
      print(ggplot2::ggplot() + ggls[[l]]$rib90 + ggls[[l]]$rib50 +
              ggls[[l]]$data + ggls[[l]]$sub + ggls[[l]]$pop +
              ggplot2::theme_bw() +
              ggplot2::theme(legend.position="none") +
              ggplot2::ggtitle(paste("Population :", l)))
    }
  }
  invisible(ggls)
}






