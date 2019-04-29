#' Create simple binned plot.

#' @param p_zstar whether to show vertical line for zstar. Default is true.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  A plot with just the frequency per bin (histogram).
#' @export



plot_hist <- function(z_vector, binv = "median", zstar,
                          binwidth, bins_l, bins_r,
                      p_title = "", p_xtitle = "z_name", p_ytitle = "Count",
                      p_maxy = NA, p_axis_txt_size = 7, p_axis_val_size = 7,
                      p_theme = "bw_light",  p_freq_color = "black",
                      p_zstar_color = "red",
                      p_freq_size = .5, p_freq_msize = 1, p_zstar_size = .5, p_zstar = T) {

    # turn data into dataframe of binned counts
    binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)

    # set theme
    if (p_theme == "bw_light") {
        p_theme <- "theme_bw() + theme_light()"
    }

    # get name of z_vector to pass as xtitle if chosen
    if (p_xtitle == "z_name") {
        p_xtitle <- deparse(substitute(z_vector))
    }
    # plot histogram
    hist_plot <- ggplot2::ggplot(binned_data, aes(y = freq_orig, x = bin)) + geom_line(colour = p_freq_color, size = p_freq_size)


    # show zstar vline if requested
    if(p_zstar) {
        hist_plot <- hist_plot +
            geom_vline(xintercept=zstar,  linetype = "solid", size = p_zstar_size, colour = p_zstar_color)
    }
    # add bpoint (so it appears in front of vline) and rest of options
    hist_plot <- hist_plot + geom_point(colour = p_freq_color, size = p_freq_msize) +
        eval(parse(text = p_theme)) +  #theme_bw() + theme_light() +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(), plot.title = element_text(hjust=0.5),
              text = element_text(size=p_axis_txt_size),
              axis.text=element_text(size=p_axis_val_size), legend.position = "none") +
        labs(title = p_title, x = p_xtitle, y = p_ytitle) +  ylim(0,p_maxy) + guides(fill=FALSE, color=FALSE)

    return(hist_plot)

    }



