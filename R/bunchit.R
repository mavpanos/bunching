#' Run the bunching estimator
#'
#' @param z_vector the vector of (unbinned) data to be analysed.
#' @param binv zstar location within its bin (min, max or median value). Default is median.
#' @param zstar the bunching point.
#' @param binwidth width of each bin.
#' @param bins_l bins to left of zstar to use.
#' @param bins_r bins to right of zstar to use.
#' @param poly order of polynomial for counterfactual fit.
#' @param bins_excl_l number of bins left of zstar to include in bunching region.
#' @param bins_excl_r number of bins right of zstar to include in bunching region.
#' @param extra_fe bin values to control for using fixed effect.
#' @param rn round number bunching (up to 2 round numbers).
#' @param n_boot number of bootstrapped iterations.
#' @param correct whether to implement correction for integration constraint.
#' @param iter_max maximum iterations for integration constraint correction.
#' @param p_title plot title.
#' @param p_xtitle plot x_axis label.
#' @param p_ytitle plot y_axis label.
#' @param p_maxy plot's maximum y_axis value.
#' @param p_txt_size text size of plot' axes' labels.
#' @param p_theme plot theme.
#' @param p_freq_color plot's frequency line color
#' @param p_cf_color plot's counterfactual line color
#' @param p_zstar_color plot's bunching region marker lines color
#' @param p_freq_size plot's frequency line thickness
#' @param p_cf_size plot's counterfactual line thickness
#' @param p_cf_msize plot's counterfactual markers' size
#' @param p_zstar_size plot's bunching region marker lines thickness
#' @param p_b Should bunching estimate be shown on plot?
#' @param p_b_xpos plot's xaxis coordinate of bunching estimate
#' @param p_b_ypos plot's yaxis coordinate of bunching estimate
#' @param p_b_size size of plot's printed bunching estimate
#' @param t0 marginal/average tax rate below zstar. see notch option.
#' @param t1 marginal/average tax rate above zstar. see notch option.
#' @param notch whether it is a notch or kink. Default is kink.
#' @param force_notch whether to enforce manual choice of zu in notch case.
#' @export

bunchit <- function(z_vector, binv = "median", zstar, binwidth, bins_l, bins_r,
                    poly, bins_excl_l, bins_excl_r, extra_fe = NA, rn = NA,
                    n_boot = 50, correct = T, iter_max = 200,
                    t0, t1, notch = F, force_notch = F,
                    p_title = "", p_xtitle = "z_name", p_ytitle = "Count",
                    p_maxy = NA, p_axis_txt_size = 7, p_axis_val_size = 7,
                    p_theme = "bw_light",  p_freq_color = "black",
                    p_cf_color = "maroon", p_zstar_color = "red",
                    p_freq_size = .5, p_cf_size = .5, p_freq_msize = 1, p_zstar_size = .5,
                    p_b = T, p_b_xpos = posx, p_b_ypos = posy, p_b_size = 3) {

    # ------------------------------------------------
    # check that inputs make sense
    # ------------------------------------------------

    # data must be a dataframe
    if(is.numeric(z_vector) == F) {
        stop("Running variable must be a numeric vector")
    }

    # binv: is it one of the 3 allowed ones?
    if(binv %in% c("min", "max", "median") == F) {
        stop("binv can only be one of 'min', 'max', 'median' ")
    }

    # zstar: is it within range? get max and min of variable
    data_varmax <- max(z_vector, na.rm = T)
    data_varmin <- min(z_vector, na.rm = T)
    if(zstar > data_varmax |  zstar < data_varmin) {
        stop("zstar (bunching point) is outside of running variable's range of values")
    }

    # binwidth: must be positive
    if(binwidth <= 0) {
        stop("Binwidth must be positive")
    }

    # bins_l: must be positive
    if(bins_l <= 0) {
        stop("bins_l must be positive")
    }

    # bins_r: must be positive
    if(bins_r <= 0) {
        stop("bins_r must be positive")
    }

    # polynomial order cannot be negative
    if(poly < 0) {
        stop("Polynomial order must be positive")
    }

    # excluded bins below zstar cannot be negative
    if(bins_excl_l < 0) {
        stop("Bins in bunching region below zstar must be a positive integer")
    }

    # excluded bins above zstar cannot be negative
    if(bins_excl_r < 0) {
        stop("Bins in bunching region above zstar must be a positive integer")
    }

    # is bunching region below zstar too wide?
    if(bins_excl_l >= bins_l - 5) {
        stop("Bunching region below zstar too wide")
    }

    # is bunching region above zstar too wide?
    if(bins_excl_r >= bins_r - 5) {
        stop("Bunching region above zstar too wide")
    }




    # check that rn does not include 0
    if( 0 %in% rn) {
        stop("Error: rn cannot include zero for a round number")
    }
    # chec that no more than two round numbers are given
    if(length(rn) > 2) {
        stop("Error: cannot have more than two unique levels for round number bunching")
    }
    # if two round numbers given, check that they are unique
    if(length(rn) == 2 & length(unique(rn)) != 2) {
        stop("Error: the two round numbers cannot be identical")
    }

    # excluded bins above zstar cannot be negative
    if(n_boot < 0) {
        stop("Bootstrap sample size must be a positive integer")
    }

    # max iteration cannot be negative
    if(iter_max < 0) {
        stop("Maximum number of iterations in integration corrrection step must be a positive integer")
    }

    # checks of inputs t0 Vs t1
    if(t0 == t1) {
        stop("Cannot calculate elasticity (t0 cannot equal t1")
    }



    # -----------------------------
    # 1. bin the data
    # -----------------------------
    # turn data into dataframe of binned counts
    binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)

    # -----------------------------
    # 2. first pass prep and fit
    # -----------------------------

    # If not a notch  or a notch but we force user choice of bins_excl_r, fit as usual

    if ((notch == F) | (notch == T & force_notch == T)) {
        # prepare data
        firstpass_prep <- bunching::prep_data_for_reg(binned_data, zstar, binwidth, bins_l, bins_r,
                                                     poly, bins_excl_l, bins_excl_r, rn, extra_fe)
        # fit firstpass model
        firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch)




        # otherwise, do data-driven notch correction
    } else if ((notch == T) & (force_notch == F)) {
        # start with only one bin above zstar
        bins_excl_r <- 1
        firstpass_prep <- bunching::prep_data_for_reg(binned_data, zstar, binwidth, bins_l, bins_r,
                                                     poly, bins_excl_l, bins_excl_r, rn, extra_fe)
        firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch)
        # extract bunching mass below and missing mass above zstar
        B_below <- firstpass$B_zl_zstar
        M_above <- -firstpass$B_zstar_zu

        # check that missing above is smaller. if not, suspicious
        if(M_above > B_below) {
            stop("Missing mass above zstar is larger than bunching mass below. Are you sure this is a notch?")
        }

        # if we reach here, we can start shifting zu bins up until B = M.
        # first, must set upper bound on number of iterations.
        #   bins above zstar cannot be more than min of:
        #       1. available bins
        #       2. remaining DoF
        available_bins <- (max(firstpass_prep$data_binned$bin) - zstar)/binwidth
        DoF_remaining <- nrow(firstpass_prep$data_binned) -  nrow(firstpass$coefficients) - 1
        notch_iterations_bound <- min(available_bins, DoF_remaining)

        # start with first bin above being bins_excl_r = 1
        zu_bin <- bins_excl_r
        while ((B_below > M_above) & (zu_bin < notch_iterations_bound)) {
            zu_bin <- zu_bin + 1
            # add next bin_excl_r dummy as column to data
            newvar <- paste0("bin_excl_r_", zu_bin)
            firstpass_prep$data_binned[[newvar]]  <- ifelse(firstpass_prep$data_binned$z_rel == zu_bin,1,0)
            # add next order bin_excl_r to formula
            firstpass_prep$model_formula <- as.formula(paste(Reduce(paste, deparse(firstpass_prep$model_formula)), newvar, sep = " + "))
            # re-fit model using the now expanded zu
            firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch)
            # get new B below and M above
            B_below <- firstpass$B_zl_zstar
            M_above <- -firstpass$B_zstar_zu
        }
        # assign final zu_bin to bins_excl_r (used for plotting)
        bins_excl_r <- zu_bin
    }

    # after fitting is done, extract info (counterfactual, residuals, etc.)
    counterfactuals_for_graph <- firstpass$cf_density
    residuals_for_boot <- firstpass$residuals
    bunchers_initial <- firstpass$bunchers_excess
    b_estimate <- firstpass$b_estimate
    e_estimate <- bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch)
    model_fit <- firstpass$coefficients


    # -----------------------------------------
    # 3. if no correction needed, do bootstrap
    # -----------------------------------------

    if(correct == F) {
        boot_results <- bunching::do_bootstrap(firstpass_prep, residuals_for_boot, boot_iterations = n_boot,
                                               correction = correct, correction_iterations = iter_max, notch = notch)
        b_sd <- boot_results$b_sd
        b_vector <- boot_results$b_vector
        e_sd <- bunching::elasticity_sd(boot_results$b_vector, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch)
        B_for_output <- bunchers_initial # this is bunchers excess. if we dont do integration constraint, this will be output
        B_sd <- boot_results$B_sd
        B_vector <- boot_results$B_vector

    }


    # ----------------------------------------------------------
    # 4. if correction needed, do that first to get residuals
    # ----------------------------------------------------------
    if (correct == T) {
        # initial correction to get vector of residuals for bootstrap for later
        firstpass_corrected <- bunching::do_correction(firstpass_prep$data_binned, firstpass_results = firstpass,
                                                       max_iterations = iter_max, notch = notch)
        b_estimate <- firstpass_corrected$b_corrected
        counterfactuals_for_graph <- firstpass_corrected$data$cf_density
        residuals_for_boot <- firstpass_corrected$data$residuals
        # we now have the correct residuals. add to our original data in firstpass_prep$data
        # applying correction each time
        boot_results <- bunching::do_bootstrap(firstpass_prep, residuals_for_boot, boot_iterations = n_boot,
                                               correction = correct, correction_iterations = iter_max, notch = notch)
        b_sd <- boot_results$b_sd
        b_vector <- boot_results$b_vector
        e_sd <- bunching::elasticity_sd(boot_results$b_vector, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch)
        B_for_output <- firstpass_corrected$B_corrected
        B_sd <- boot_results$B_sd
        B_vector <- boot_results$B_vector
        model_fit <- firstpass_corrected$coefficients
    }


    # ---------------------------------------------
    # 5. make plot
    # ---------------------------------------------
    # get max of binned_data to position b
    zmax <- max(firstpass_prep$data_binned$bin)
    posx <- zstar + (zmax - zstar)*.7
    posy <- max(firstpass_prep$data_binned$freq_orig, counterfactuals_for_graph)*.8

    # get name of z_vector to pass as xtitle if chosen
    if (p_xtitle == "z_name") {
        p_xtitle <- deparse(substitute(z_vector))
    }

    # theme
    if (p_theme == "bw_light") {
        p_theme <- "theme_bw() + theme_light()"
    }
    p <- bunching::plot_bunching(firstpass_prep$data_binned, cf = counterfactuals_for_graph, zstar,
                   binwidth, bins_excl_l, bins_excl_r,
                   p_title, p_xtitle, p_ytitle, p_maxy, p_axis_txt_size, p_axis_val_size,
                   p_theme, p_freq_color, p_cf_color, p_zstar_color,
                   p_freq_size, p_cf_size, p_freq_msize, p_zstar_size,
                   p_b, b = b_estimate, b_sd = b_sd,
                   p_b_xpos, p_b_ypos, p_b_size)


    output <- list("data" = firstpass_prep$data_binned,
                   "cf" = counterfactuals_for_graph,
                   "B" = B_for_output,
                   "B_vector" = B_vector,
                   "B_sd" = B_sd,
                   "b" = b_estimate,
                   "b_vector" = b_vector,
                   "b_sd" = b_sd,
                   "e" = e_estimate,
                   "e_sd" = e_sd,
                   "plot" = p,
                   "modelfit" = model_fit)
    return(output)
}

