#' Bunching Plot
#'
#' Creates the bunching plot.

#'
#' @param binned_data binned data with frequency and estimated counterfactual.
#' @param b normalized bunching estimate.
#' @param b_sd standard deviation of the normalized bunching estimate.
#' @param e elasticity estimate.
#' @param e_sd standard deviation of the elasticity estimate.
#' @param cf the counterfactual to be plotted.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  \code{plot_bunching} returns a plot with the frequency, counterfactual and bunching region demarcated. Can also include the bunching and elasticity estimate if specified.
#'
#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' prepped_data <- prep_data_for_fit(binned_data, zstar = 10000, binwidth = 50,
#'                                   bins_l = 20, bins_r = 20, poly = 4)
#' fitted <- fit_bunching(thedata = prepped_data$data_binned,
#'                        themodelformula = prepped_data$model_formula)
#' plot_bunching(z_vector = bunching_data$kink_vector,
#'               binned_data = prepped_data$data_binned,
#'               cf = fitted$cf_density, zstar = 10000,
#'               binwidth = 50, bins_excl_l = 0 , bins_excl_r = 0,
#'               b = 1.989, b_sd = 0.005, p_b = TRUE)
#' @export
#'
plot_bunching <- function(z_vector, binned_data, cf, zstar,
                          binwidth, bins_excl_l = 0, bins_excl_r = 0,
                          p_title = "", p_xtitle = deparse(substitute(z_vector)), p_ytitle = "Count",
                          p_miny = 0, p_maxy = NA, p_ybreaks = NA,
                          p_title_size = 11, p_axis_title_size = 10, p_axis_val_size = 8.5,
                          p_freq_color = "black", p_cf_color = "maroon",
                          p_zstar_color = "red", p_grid_major_y_color = "lightgrey",
                          p_freq_size = .5, p_freq_msize = 1, p_cf_size = .5, p_zstar_size = .5,
                          p_b = FALSE, b = NA, b_sd = NA, p_e = FALSE, e = NA, e_sd = NA, p_b_e_xpos = NA,
                          p_b_e_ypos = NA, p_b_e_size = 3,
                          t0 = NA, t1 = NA, notch = FALSE,
                          p_domregion_color = NA, p_domregion_ltype = NA) {


    # if p_b_e_xpos/p_y_e_xpos not chosen, set them. get data ranges
    zmin <- min(binned_data$bin, na.rm = TRUE)
    zmax <- max(binned_data$bin,na.rm = TRUE)     # also needed to customize plot region
    maxy <- max(binned_data$freq_orig, cf)

    # bunching/elasticity estimates' x position
    if(is.na(p_b_e_xpos)) {
        if (notch == TRUE) {
            p_b_e_xpos <- zmin + (zstar - zmin)*.3
        } else {
            p_b_e_xpos <- zstar + (zmax - zstar)*.7
        }
    }

    # y position
    if(is.na(p_b_e_ypos)) {
        p_b_e_ypos <- maxy * .8
    }

    # prepare data for graphing in "long" structure
    df_plot <- data.frame(cbind(binned_data[,c("bin", "freq_orig")], "cf_graph" = cf))
    df_plot <- df_plot %>%  tidyr::gather(key, value, freq_orig, cf_graph)

    # get position of bunching region lower and upper bounds and add them to vector that includes kinklevel
    lb <- zstar - bins_excl_l * binwidth
    ub <- zstar + bins_excl_r * binwidth
    # create vector of linetypes (solid Vs dashed) for each
    vlines <- c(lb, zstar, ub)
    vlines_type <- rep("solid", 3) # make them all solid. if one is different, make dashed
    vlines_type[vlines != zstar] <- "dashed"
    vlines_color <- rep(p_zstar_color, 3)

    # if it's a notch, add domregion bin marker to vlines
    if(notch == TRUE) {
        bin_domregion <- domregion(zstar, t0, t1, binwidth)$zD
        vlines <- c(vlines, bin_domregion)
        vlines_type <- c(vlines_type, p_domregion_ltype)
        vlines_color <- c(vlines_color, p_domregion_color)
    }

    # combine b and b_sd (e and e_sd) to pass to graph (if b_sd is not NA, otherwise only report point estimate)
    if(!is.na(b_sd)) {
        b_estimates <- paste0("b = ", round(b,3), "(", round(b_sd,3),  ")")
        e_estimates <- paste0("e = ", round(e,3), "(", round(e_sd, 3),  ")")
    } else {
        b_estimates <- paste0("b = ", round(b,3))
        e_estimates <- paste0("e = ", round(e,3))
    }

    # prepare color vector for main lines
    freq_cf_colors <- c("freq_orig" = p_freq_color, "cf_graph" = p_cf_color)
    freq_cf_sizes <- c("freq_orig" = p_freq_size, "cf_graph" = p_cf_size)

    # the plot
    bunch_plot <- ggplot2::ggplot(df_plot, aes(y = value, x = bin, color = key, size = key)) +
        # plot main lines, pass size and color options
        ggplot2::geom_line() + ggplot2::scale_size_manual(values = freq_cf_sizes) + ggplot2::scale_colour_manual(values = freq_cf_colors) +
        # plot vertical lines first so that they appear behind plot binpoints
        ggplot2::geom_vline(xintercept=vlines,  linetype=vlines_type,  color = vlines_color, size=p_zstar_size) +
        # plot connector point marker for freq line
        ggplot2::geom_point(data = df_plot[df_plot$key == "freq_orig",], size = p_freq_msize) +
        theme_classic() + #eval(parse(text = p_theme))
        ggplot2::theme(panel.grid.major.y = element_line(colour = p_grid_major_y_color),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       plot.title = element_text(hjust = 0.5, size = p_title_size),
                       text = element_text(size = p_axis_title_size),
                       axis.text=element_text(size = p_axis_val_size),
                       legend.position = "none",
                       panel.background = element_blank(),
                       panel.grid = element_blank(),
                       panel.border = element_blank()) +
        ggplot2::labs(title = p_title, x = p_xtitle, y = p_ytitle) +
        ggplot2::guides(fill = FALSE, color = FALSE)

    # pass choice of ylim and p_ybreaks
    if(sum(is.na(p_ybreaks)) == length(p_ybreaks)) {
        bunch_plot <-  bunch_plot + scale_y_continuous(limits= c(p_miny, p_maxy))
    } else {
        bunch_plot <- bunch_plot + scale_y_continuous(limits= c(p_miny, p_maxy), breaks = p_ybreaks)
    }

    # choice to show b (and e) on plot or not
    if(p_b == TRUE) {
        if(p_e == TRUE) {
            # if both TRUE, add both on separate lines
            text_to_print <- paste0(b_estimates, "\n", e_estimates)
        } else { # add only b
            text_to_print <- b_estimates
        }
        bunch_plot <- bunch_plot +
            ggplot2::annotate("text", x = p_b_e_xpos, y = p_b_e_ypos,
                              label = text_to_print, size = p_b_e_size)
    }

    return(bunch_plot)
}

