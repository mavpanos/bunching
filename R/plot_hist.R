#'Plot Histogram
#'
#' Create a binned plot for quick exploration without estimating bunching mass.

#' @param p_zstar whether to show vertical line for zstar. Default is TRUE.
#' @inheritParams bunchit
#' @return  \code{plot_hist} returns a list with the following:
#'   \item{plot}{the plot of the density without estimating a counterfactual.}
#'   \item{data}{the binned data used for the plot.}

#' @seealso \code{\link{bunchit}}
#'
#' @examples
#'
#' # visualize a distribution
#' data(bunching_data)
#' plot_hist(z_vector = bunching_data$kink_vector,
#' binv = "median", zstar = 10000,
#' binwidth = 50, bins_l = 40, bins_r = 40)$plot
#'
#' @export




plot_hist <- function(z_vector, binv = "median", zstar, binwidth, bins_l, bins_r,
                      p_title = "", p_xtitle = "z_name", p_ytitle = "Count",
                      p_title_size = 11, p_axis_title_size = 10, p_axis_val_size = 8.5,
                      p_miny = 0, p_maxy = NA, p_ybreaks = NA,
                      p_grid_major_y_color = "lightgrey",
                      p_freq_color = "black", p_zstar_color = "red",
                      p_freq_size = .5, p_freq_msize = 1, p_zstar_size = .5, p_zstar = TRUE) {

    # ------------------------------------------------
    #               0. check inputs
    # ------------------------------------------------

    # z_vector must be numeric
    if(!is.numeric(z_vector)) {
        stop("z_vector must be a numeric vector")
    }

    # binv: is it one of the 3 allowed ones?
    if(binv %in% c("min", "max", "median") == FALSE) {
        stop("binv can only be one of 'min', 'max', 'median'")
    }

    # zstar: is it within range? get max and min of variable
    data_varmax <- max(z_vector, na.rm = TRUE)
    data_varmin <- min(z_vector, na.rm = TRUE)
    if(zstar > data_varmax |  zstar < data_varmin) {
        stop("zstar is outside of z_vector's range of values")
    }

    # binwidth: must be positive, numeric
    if(binwidth <= 0 | !is.numeric(binwidth)) {
        stop("Binwidth must be a positive number")
    }

    # check that bins_l are positive and integers
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

    # bins_l: must be positive
    if(bins_l <= 0 | !is.wholenumber(bins_l)) {
        stop("bins_l must be a positive integer")
    }

    # bins_r: must be positive
    if(bins_r <= 0 | !is.wholenumber(bins_r)) {
        stop("bins_r must be a positive integer")
    }


    # is p_title a string?
    if(!is.character(p_title)) {
        stop("p_title must be a string")
    }

    # is p_xtitle a string?
    if(!is.character(p_xtitle)) {
        stop("p_xtitle must be a string")
    }

    # is p_ytitle a string?
    if(!is.character(p_ytitle)) {
        stop("p_ytitle must be a string")
    }

    # is p_title_size numeric and positive?
    if(p_title_size <= 0 | !is.numeric(p_title_size)) {
        stop("p_title_size must be a positive numeric value")
    }

    # is p_axis_title_size numeric and positive?
    if(p_axis_title_size <= 0 | !is.numeric(p_axis_title_size)) {
        stop("p_axis_title_size must be a positive numeric value")
    }

    # is p_axis_val_size numeric and positive?
    if(p_axis_val_size <= 0 | !is.numeric(p_axis_val_size)) {
        stop("p_axis_val_size must be a positive numeric value")
    }

    # is p_miny numeric?
    if(!is.na(p_miny) & !is.numeric(p_miny)) {
        stop("plot's minimum y-value p_miny must be numeric")
    }

    # is p_maxy numeric?
    if(!is.na(p_maxy) & !is.numeric(p_maxy)) {
        stop("plot's maximum y-value p_miny must be numeric")
    }

    # is p_maxy > p_miny?
    if(!is.na(p_maxy) &  (p_maxy < p_miny)) {
        stop("p_maxy cannot be smaller than p_miny")
    }

    # are all elements in p_ybreaks numeric (if specified, otherwise NA)
    # if more than one entered, check that all are numeric
    if(length(p_ybreaks) > 1) {
        # if at least one has NA, stop
        if(length(p_ybreaks[is.character(p_ybreaks) | is.na(p_ybreaks)]) > 0) {
            stop("p_ybreaks must only contain numeric values")
        }
    } else {
        # if only of length 1, and not NA nor numeric,
        if(!is.na(p_ybreaks) & !is.numeric(p_ybreaks)) {
            stop("p_ybreaks must only contain numeric values")
        }
    }

    # is p_grid_major_y_color a string?
    if(!is.character(p_grid_major_y_color)) {
        stop("p_grid_major_y_color must a string, e.g. 'blue'")
    }

    # is p_freq_color a string?
    if(!is.character(p_freq_color)) {
        stop("p_freq_color must be a string, e.g. 'black'")
    }

    # is p_zstar_color a string?
    if(!is.character(p_zstar_color)) {
        stop("p_zstar_color must be a string, e.g. 'red'")
    }

    # is p_freq_size numeric?
    if(!is.numeric(p_freq_size)) {
        stop("p_freq_size must be numeric")
    }

    # is p_freq_msize numeric and positive?
    if(p_freq_msize <= 0 | !is.numeric(p_freq_msize)) {
        stop("p_freq_msize must be a positive numeric value")
    }

    # is p_zstar logical?
    if(!is.logical(p_zstar)) {
        stop("p_zstar can only be TRUE or FALSE")
    }

    # is p_zstar_size numeric?
    if(p_zstar_size <= 0 | !is.numeric(p_zstar_size)) {
        stop("p_zstar_size must be a positive numeric value")
    }


    # ------------------------------------------------
    #               1. bin the data
    # ------------------------------------------------
    # convert z_vector into dataframe of binned counts
    binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)
    #drop extra freq and z columns (uneccesary for hist)
    binned_data$freq <- NULL
    binned_data$z <- NULL

    # ------------------------------------------------
    #               2. make plot
    # ------------------------------------------------


    # get name of z_vector to pass as xtitle if chosen
    if (p_xtitle == "z_name") {
        p_xtitle <- deparse(substitute(z_vector))
    }
    # plot basic histogram
    hist_plot <- ggplot2::ggplot(binned_data, aes(y = freq_orig, x = bin)) + geom_line(colour = p_freq_color, size = p_freq_size)

    # show zstar vline if requested through p_zstar
    if(p_zstar) {
        hist_plot <- hist_plot +
            ggplot2::geom_vline(xintercept=zstar,  linetype = "solid", size = p_zstar_size, colour = p_zstar_color)
    }
    # add zstar (so it appears in front of vline) and rest of options
    hist_plot <- hist_plot + ggplot2::geom_point(colour = p_freq_color, size = p_freq_msize) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(colour = p_grid_major_y_color),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size = p_title_size),
              text = element_text(size = p_axis_title_size),
              axis.text = element_text(size = p_axis_val_size),
              legend.position = "none",
              panel.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank()) +
        ggplot2::labs(title = p_title, x = p_xtitle, y = p_ytitle) +
        ggplot2::guides(fill = FALSE, color = FALSE)

    # pass choice of ylim and p_ybreaks
    if(sum(is.na(p_ybreaks)) == length(p_ybreaks)) {
        hist_plot <-  hist_plot + scale_y_continuous(limits= c(p_miny, p_maxy))
    } else {
        hist_plot <- hist_plot + scale_y_continuous(limits= c(p_miny, p_maxy), breaks = p_ybreaks)
    }

    return(list(data = binned_data,
                plot = hist_plot))
}



