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
plot_spline_v2 <- function(model, limits = NULL, plot_which = NULL, plot_type = "mean",
                        shade = FALSE, silent = FALSE) {


    ## Rename the input data to prevent future accidental name change.
    data <- model$data
    colnames(data) <- c("x", "y", "sub", "pop")

    lvl_pop <- levels(data$pop)
    subs_in_pop <- tapply(as.character(data$sub), data$pop, unique, simplify = FALSE)

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

    ## if (shade) {
    ##     if (plot_which != "dme__") {
    ##         stop("shade only works for models with single population.")
    ##     }
    ##     ## thin the sample by 10
    ##     thin_idx <- seq(10, NCOL(model$samples$population), 10)
    ##     thin_y <- get_plot_y(model$samples$population[, thin_idx])
    ##     colnames(thin_y) <- thin_idx
    ##     prob <- c(0.25, 0.75)
    ##     prob_outer <- c(0.05, 0.95)
    ##     plotdat_thin <- apply(thin_y, 1, function(x) quantile(x, c(prob, prob_outer))) %>%
    ##         {tibble::as_tibble(t(.))} %>%
    ##         dplyr::mutate(x = plot_x)
    ## }

    ## Reformat dataframe for ggplot

    EPS <- 1e-6
    fine <- 200                         # how fine the plot_x should be?

    plot_x <- c(seq(min(limits), max(limits), length.out = fine),
                model$pop$basis$knots, model$sub$basis$knots) %>%
        unique() %>% sort()

    get_model_mat <- function(plot_x, binfo) {

        knots <- binfo$knots
        names(knots) <- NULL
        deg <- binfo$degree

        if (binfo$type == "tpf") {
            model_mat <- get_design_tpf(plot_x, knots, deg)$design
        } else if (binfo$type == "bs-ridge") {
            model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                               outer.ok = TRUE)
        } else if (binfo$type == "bs") {
            model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                               outer.ok = TRUE)
            model_mat <- model_mat %*% binfo$trans_mat
        } else {
            stop("Unknown type of basis.")
        }
        model_mat
    }

    model_mat_pop <- get_model_mat(plot_x, model$basis$pop)
    model_mat_sub <- get_model_mat(plot_x, model$basis$sub)
    if (NROW(model_mat_pop) != NROW(model_mat_sub)) {
        stop("Dimension error.")
    }

    ggls <- list()
    for (l in plot_which) {
        subs_in_l <- subs_in_pop[[l]]
        plot_y_pop <- as.numeric(model_mat_pop %*% coef_pop[, l])
        plot_y_sub <- plot_y_pop + (model_mat_sub %*% dev_sub[, subs_in_l])
        plotdat_pop <- tibble::tibble(x = plot_x, y = plot_y_pop)
        plotdat_sub <- tibble::as_tibble(plot_y_sub) %>%
            dplyr::mutate(x = plot_x) %>%
            tidyr::gather("sub", "y", -"x")
        ggls$data[[l]] <- geom_point(aes_(~x, ~y, col = ~sub), data[data$pop == l, ])
        ## if (shade) {
        ##     ggls$rib90 <- geom_ribbon(aes_(~x, ymin = ~`5%`, ymax = ~`95%`),
        ##                                 plotdat_thin, fill = "grey85", alpha = 0.7)
        ##     ggls$rib50 <- geom_ribbon(aes_(~x, ymin = ~`25%`, ymax = ~`75%`),
        ##                               plotdat_thin, fill = "grey65", alpha = 0.5)
        ## }
        ggls$sub[[l]] <- geom_line(aes_(~x, ~y, col = ~sub, group = ~sub), plotdat_sub)
        ggls$pop[[l]] <- geom_line(aes_(~x, ~y), data = plotdat_pop)

        if (!silent) {
            print(ggplot2::ggplot() + ggls$data[[l]] + ggls$pop[[l]] + ggls$sub[[l]] +
                  ggplot2::theme_bw() +
                  ggplot2::theme(legend.position="none"))
        }
    }
    ggls
}








