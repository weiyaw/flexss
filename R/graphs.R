################################################################
##########################            ##########################
##########################    PLOT    ##########################
##########################            ##########################
################################################################

## Plot the population and subject curves (suitable for models with multiple populations)
## 'plot_which' takes a vector of names (for models with multiple populations)
## type can be either of "mean", "mle" or "map"
## If shade is true, thin the population samples by 10 and plot the thinned samples.

## Requires: model$basis$knots (including extrema), model$basis$type,
##   model$basis$degree, model$info$lvl_pop, model$data, model$mean,
##   model$samples$population (if shade is true)

#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes_ geom_point geom_line
plot_spline <- function(model, limits = NULL, plot_which = NULL, plot_type = "mean",
                        shade = FALSE, silent = FALSE) {

    EPS <- 1e-6
    fine <- 200                         # how fine the plot_x should be?
    knots <- model$basis$knots
    names(knots) <- NULL
    type <- model$basis$type
    deg <- model$basis$degree

    ## Rename the input data to prevent future accidental name change.
    data <- model$data
    colnames(data) <- c("x", "y", "sub", "pop")

    lvl_pop <- levels(data$pop)

    ## Use the range of the predictor if 'limits' is missing
    if (is.null(limits)) {
        limits <- range(data$x)
    } else if (!(is.vector(limits) && length(limits) == 2)) {
        stop("'limits' must be a vector of length 2.")
    }

    ## Plot all population if 'plot_which' is missing
    if (is.null(plot_which)) {
        plot_which <- lvl_pop
    } else if (!is.vector(plot_which)) {
        stop("'plot_which' must be a vector.")
    } else if (!all(plot_which %in% lvl_pop)) {
        stop("Invalid population levels in 'plot_which'.")
    }

    ## 'x' axis of the plot
    knots_within <- knots[knots > (min(limits) - EPS) &
                          knots < (max(limits) + EPS)]
    plot_x <- unique(c(seq(min(limits), max(limits), length.out = fine),
                       knots_within))
    plot_x <- plot_x[order(plot_x)]

    ## Model matrix
    if (type == "tpf") {
        model_mat <- get_design_tpf(plot_x, knots, deg)$design
    } else if (type == "bs-ridge") {
        model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                           outer.ok = TRUE)
    } else if (type == "bs") {
        model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                           outer.ok = TRUE)
        model_mat <- model_mat %*% model$basis$trans_mat
    } else {
        stop("Unknown type of model.")
    }

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

    ## Helper function to calculate the y axis of the plot
    get_plot_y <- function(coef) {
        model_mat %*% coef
    }

    ## y axis of the plot
    if (length(lvl_pop) == 0) {
        coef_sub <- dev_sub + coef_pop
        ## Align the data format of a single population model to a multiple
        ## population model. Create a dummy name for the population
        plot_which <- "dummy"
        data$pop <- "dummy"
        plot_y_pop <- list(dummy = get_plot_y(coef_pop))
        plot_y_sub <- list(dummy = get_plot_y(coef_sub))
    } else {
        coef_sub <- mapply(`+`, dev_sub[plot_which], coef_pop[plot_which],
                           SIMPLIFY = FALSE)
        plot_y_pop <- lapply(coef_pop[plot_which], get_plot_y)
        plot_y_sub <- lapply(coef_sub[plot_which], get_plot_y)
    }

    if (shade) {
        if (plot_which != "dummy") {
            stop("shade only works for models with single population.")
        }
        ## thin the sample by 10
        thin_idx <- seq(10, NCOL(model$samples$population), 10)
        thin_y <- get_plot_y(model$samples$population[, thin_idx])
        colnames(thin_y) <- thin_idx
        prob <- c(0.25, 0.75)
        prob_outer <- c(0.05, 0.95)
        plotdat_thin <- apply(thin_y, 1, function(x) quantile(x, c(prob, prob_outer))) %>%
            {tibble::as_tibble(t(.))} %>%
            dplyr::mutate(x = plot_x)
    }
    ## Reformat dataframe for ggplot
    plotdat_pop <- list()
    plotdat_sub <- list()
    for (i in plot_which) {
        plotdat_pop[[i]] <- data_frame(x = plot_x, y = as.numeric(plot_y_pop[[i]]))
        plotdat_sub[[i]] <- tibble::as_tibble(plot_y_sub[[i]]) %>%
            dplyr::mutate(x = plot_x) %>%
            tidyr::gather("sub", "y", -"x")
         }

    ggls <- list()
    for (i in plot_which) {
        ggls$data <- geom_point(aes_(~x, ~y, col = ~sub), data[data$pop == i, ])
        if (shade) {
            ggls$rib90 <- geom_ribbon(aes_(~x, ymin = ~`5%`, ymax = ~`95%`),
                                        plotdat_thin, fill = "grey85", alpha = 0.7)
            ggls$rib50 <- geom_ribbon(aes_(~x, ymin = ~`25%`, ymax = ~`75%`),
                                      plotdat_thin, fill = "grey65", alpha = 0.5)
        }
        ggls$sub <- geom_line(aes_(~x, ~y, col = ~sub, group = ~sub), plotdat_sub[[i]])
        ggls$pop <- geom_line(aes_(~x, ~y), data = plotdat_pop[[i]])

        if (!silent) {
            print(ggplot2::ggplot() + ggls + theme_bw() + theme(legend.position="none"))
        }
    }
    ggls
}








