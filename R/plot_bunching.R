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
#' @export


plot_bunching <- function(binned_data, cf, zstar,
                          binwidth, bins_excl_l, bins_excl_r,
                          p_title, p_xtitle, p_ytitle, p_miny, p_maxy, p_ybreaks,
                          p_title_size, p_axis_title_size, p_axis_val_size,
                          p_freq_color, p_cf_color, p_zstar_color, p_grid_major_y_color,
                          p_freq_size, p_cf_size, p_freq_msize, p_zstar_size,
                          p_b, b, b_sd, p_e, e, e_sd, p_b_e_xpos, p_b_e_ypos, p_b_e_size,
                          t0 = NA, t1 = NA, notch = F,
                          p_domregion_color = NA, p_domregion_ltype = NA, n_boot) {

    # get upper bound to customize plot region
    zmax <- max(binned_data$bin,na.rm = T)

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
    if(notch == T) {
        bin_domregion <- domregion(zstar, t0, t1, binwidth)$zD
        vlines <- c(vlines, bin_domregion)
        vlines_type <- c(vlines_type, p_domregion_ltype)
        vlines_color <- c(vlines_color, p_domregion_color)
    }

    # combine b and b_sd (e and e_sd) to pass to graph (if n_boot == 0, only report point estimate)
    if(n_boot > 0) {
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
                       plot.title = element_text(hjust=0.5, size = p_title_size),
                       text = element_text(size=p_axis_title_size),
                       axis.text=element_text(size=p_axis_val_size),
                       legend.position = "none",
                       panel.background = element_blank(),
                       panel.grid=element_blank(),
                       panel.border=element_blank()) +
        ggplot2::labs(title = p_title, x = p_xtitle, y = p_ytitle) +
        ggplot2::guides(fill=FALSE, color=FALSE)

    # pass choice of ylim and p_ybreaks
    if(sum(is.na(p_ybreaks)) == length(p_ybreaks)) {
        bunch_plot <-  bunch_plot + scale_y_continuous(limits=c(p_miny, p_maxy))
    } else {
        bunch_plot <- bunch_plot + scale_y_continuous(limits=c(p_miny, p_maxy), breaks = p_ybreaks)
    }

    # choice to show b (and e) on plot or not
    if(p_b == T) {
        if(p_e == T) {
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

