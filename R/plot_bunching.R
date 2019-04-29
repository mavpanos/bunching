#' Create bunching plot.

#'
#' @param binned_data binned data with frequency and counterfactual estimated
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  A plot with the frequency, counterfactual and bunching region demarcated. Can also include the bunching estimate in the plot if specified.
#' @export


plot_bunching <- function(binned_data, cf, zstar,
                          binwidth, bandwidth, bins_excl_l, bins_excl_r,
                          p_title, p_xtitle, p_ytitle, p_maxy, p_txt_size,
                          p_theme, p_freq_color, p_cf_color, p_zstar_color,
                          p_freq_size, p_cf_size, p_cf_msize, p_zstar_size,
                          p_b, b, b_sd,
                          p_b_xpos, p_b_ypos, p_b_size) {

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

    # combine b and b_sd to pass to graph
    b_estimates <- paste0("b = ", sprintf("%.3f", b), "(", b_sd, ")")

    # the plot
    bunch_plot <- ggplot2::ggplot(df_plot, aes(y = value, x = bin, color = key)) +
        # plot vertical lines first so that they appear behind plot binpoints
        geom_vline(xintercept=vlines,  linetype=vlines_type,  color = p_zstar_color, size=p_zstar_size) +
        # plot main graph lines
        geom_point(data = df_plot[df_plot$key == "freq_orig",], size = p_cf_msize) + geom_line(size = p_cf_size) +
        geom_line(data = df_plot[df_plot$key == "cf_graph",], size = p_freq_size) +
        scale_colour_manual(values=c(p_cf_color, p_freq_color)) +
        eval(parse(text = p_theme)) +  #theme_bw() + theme_light() +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(), plot.title = element_text(hjust=0.5),
              text = element_text(size=p_txt_size)) +
        labs(title = p_title, x = p_xtitle, y = p_ytitle) +  ylim(0,p_maxy) + guides(fill=FALSE, color=FALSE)

    # choice to show b on plot or not
    if(p_b == T) {
        bunch_plot <- bunch_plot +
            annotate("text", x = p_b_xpos, y = p_b_ypos,
                     label = b_estimates, size = p_b_size)
    }

    return(bunch_plot)
}
