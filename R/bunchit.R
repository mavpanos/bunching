#' Bunching Estimator
#'
#' Implement the bunching estimator in a kink or notch setting.
#'
#' @param z_vector a numeric vector of (unbinned) data to be analysed.
#' @param binv a string setting location of zstar within its bin ("min", "max" or "median" value). Default is median.
#' @param zstar a numeric value for the the bunching point.
#' @param binwidth a numeric value for the width of each bin.
#' @param bins_l number of bins to left of zstar to use in analysis.
#' @param bins_r number of bins to right of zstar to use in analysis.
#' @param poly a numeric value for the order of polynomial for counterfactual fit.
#' @param bins_excl_l number of bins to left of zstar to include in bunching region.
#' @param bins_excl_r number of bins to right of zstar to include in bunching region.
#' @param extra_fe a numeric vector of bin values to control for using fixed effects.
#' @param rn a numeric vector of round numbers (up to 2) to control for.
#' @param n_boot number of bootstrapped iterations.
#' @param correct if TRUE, implements correction for integration constraint.
#' @param iter_max maximum iterations for integration constraint correction.
#' @param p_title plot title.
#' @param p_xtitle plot x_axis label.
#' @param p_ytitle plot y_axis label.
#' @param p_miny plot's minimum y_axis value. Default is zero.
#' @param p_maxy plot's maximum y_axis value.
#' @param p_ybreaks plot's sequence of values to use for y-axis grid.
#' @param p_axis_title_size size of plot's axes' title labels.
#' @param p_axis_val_size size of plot's axes' numeric labels.
#' @param p_theme plot theme (in ggplot2 format).
#' @param p_freq_color plot's frequency line color.
#' @param p_cf_color plot's counterfactual line color.
#' @param p_zstar_color plot's bunching region marker lines color.
#' @param p_grid_major_y_color plot's y-axis major grid line color.
#' @param p_freq_size plot's frequency line thickness.
#' @param p_cf_size plot's counterfactual line thickness.
#' @param p_freq_msize plot's frequency line marker size.
#' @param p_zstar_size plot's bunching region marker lines thickness.
#' @param p_b if TRUE, plot also includes bunching estimate. Default is TRUE.
#' @param p_e if TRUE, plot also includes elasticity estimate. Only shown if bunching estimate shown. Default is TRUE.
#' @param p_b_xpos plot's xaxis coordinate of bunching estimate.
#' @param p_b_ypos plot's yaxis coordinate of bunching estimate.
#' @param p_b_size size of plot's printed bunching estimate.
#' @param t0 numeric value between 0 and 1 setting the marginal/average tax rate below zstar, depending on kink/notch choice. see notch parameter.
#' @param t1 numeric value between 0 and 1 setting the marginal/average tax rate above zstar, depending on kink/notch choice. see notch parameter.
#' @param notch if TRUE, zstar treated as notch. Default is kink.
#' @param force_notch if TRUE, user choice of zu (upper limit of bunching region) is enforced. Default is FALSE (zu set by setting bunching equal to missing mass).
#' @param e_parametric if TRUE, elasticity is estimated using parametric specification (quasi-linear and iso-elastic utility function). Default is FALSE (which estimates reduced-form approximation).
#' @param p_domregion_color plot's dominated region marker line color (notch case).
#' @param p_domregion_ltype a string for the vertical line type marking the dominated region (zD) in the plot (notch case only).
#' @param seed a numeric value for bootstrap seed (random re-sampling of residuals).

#' @details bunchit implements the bunching estimator in both kink and notch settings. It bins the running variable, fits a counterfactual density, and estimates the bunching mass (normalized and not), the elasticity and the location of the marginal buncher. In the case of notches, it also finds the dominated region and estimates the fraction of observations located in it.

#' @return \code{bunchit} returns a list of results, both for visualizing and for further analysis of the data underlying the estimates. These include:
#'   \item{plot}{The bunching plot.}
#'   \item{data}{The binned data used for estimation.}
#'   \item{cf}{The estimated counterfactuals.}
#'   \item{B}{The estimated excess mass (not normalized).}
#'   \item{B_vector}{The vector of bootstrapped B's.}
#'   \item{B_sd}{The standard deviation of B_vector.}
#'   \item{b}{The estimated excess mass (normalized).}
#'   \item{b_vector}{The vector of bootstrapped b's.}
#'   \item{b_sd}{The standard deviation of b_vector.}
#'   \item{e}{The estimated elasticity.}
#'   \item{e_vector}{The vector of bootstrapped elasticities (e).}
#'   \item{e_sd}{The standard deviation of e_vector.}
#'   \item{alpha}{The estimated fraction of bunchers in dominated region (notch case).}
#'   \item{alpha_vector}{The vector of bootstrapped alphas.}
#'   \item{alpha_sd}{The standard deviation of alpha_vector.}
#'   \item{model_fit}{The model fit on the actual (i.e. not bootstrapped) data.}
#'   \item{zD}{The value demarcating the dominated region (notch case).}
#'   \item{zD_bin}{The bin above zstar demarcating the dominated region (notch case).}
#'   \item{marginal_buncher}{The location (z value) of the marginal buncher.}
#'   \item{marginal_buncher_vector}{The vector of bootstrapped marginal_buncher values.}
#'   \item{marginal_buncher_sd}{The standard deviation of marginal_buncher_vector.}
#' @import BB
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @seealso \code{\link{plot_hist}}
#' @export

bunchit <- function(z_vector, binv = "median", zstar, binwidth, bins_l, bins_r,
                    poly, bins_excl_l, bins_excl_r, extra_fe = NA, rn = NA,
                    n_boot = 50, correct = TRUE, iter_max = 200,
                    t0, t1, notch = FALSE, force_notch = FALSE, e_parametric = FALSE,
                    p_title = "", p_xtitle = "z_name", p_ytitle = "Count",
                    p_axis_title_size = 7, p_axis_val_size = 7, p_miny = 0, p_maxy = NA, p_ybreaks = NULL,
                    p_theme = "theme_classic()",  p_freq_color = "black",
                    p_cf_color = "maroon", p_zstar_color = "red", p_grid_major_y_color = "lightgrey",
                    p_freq_size = .5, p_freq_msize = 1, p_cf_size = .5, p_zstar_size = .5,
                    p_b = TRUE, p_e = TRUE, p_b_xpos = NA, p_b_ypos = NA, p_b_size = 3,
                    p_domregion_color = "blue", seed = NA, p_domregion_ltype="longdash") {

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

    # zstar cannot be zero (elasticity needs Dz/zstar estimate)
    if(zstar == 0) {
        stop("zstar (bunching point) cannot be zero. If this your true bunching point, must re-centre it away from zero.")
    }

    # binwidth: must be positive
    if(binwidth <= 0) {
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

    # polynomial order cannot be negative
    if(poly < 0 | !is.wholenumber(poly)) {
        stop("Polynomial order must be a non-negative integer")
    }

    # excluded bins below zstar cannot be negative
    if(bins_excl_l < 0 | !is.wholenumber(bins_excl_l)) {
        stop("Number of bins in bunching region below zstar must be a non-negative integer")
    }

    # excluded bins above zstar cannot be negative
    if(bins_excl_r < 0 | !is.wholenumber(bins_excl_r)) {
        stop("Number of bins in bunching region above zstar must be a non-negative integer")
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
        stop("Error: rn cannot include zero as a round number")
    }

    # if not NA, are all round numbers integers?
    if(sum(is.na(rn)) == 0 & sum(is.wholenumber(rn)) != length(rn)) {
        stop("Round number(s) must be integer(s)")
    }

    # check that no more than two round numbers are given
    if(length(rn) > 2) {
        stop("Error: cannot have more than two unique levels for round number bunching")
    }

    # if two round numbers given, check that they are unique
    if(length(rn) == 2 & length(unique(rn)) != 2) {
        stop("Error: the two round numbers cannot be identical")
    }

    # if rn is not NA, are round numbers included in data range?
    if(sum(is.na(rn)) == 0 & sum(rn > data_varmax) != 0) {
        stop("rn includes round numbers outside of running variable's range of values")
    }

    # if rn is not NA,
    # is data range too small to allow for round numbers? (e.g. data range is 0-400 but rn = 500)
    if(sum(is.na(rn)) == 0 & sum(rn > (data_varmax-data_varmin))) {
        stop("rn includes round numbers that are too large for your running variable's range of values. \n Use smaller values")
    }

    # number of bootstrap iterations must be non-negative integer
    if(n_boot < 0 | !is.wholenumber(n_boot)) {
        stop("Bootstrap sample size must be a non-negative integer")
    }

    # flag if bootstrap samples less than 100
    if(n_boot > 0 & n_boot < 100) {
        warning(paste0("You chose n_boot = ", n_boot, ". Are you sure this is large enough?"))
    }


    # correct must be TRUE/FALSE
    if(!is.logical(correct)) {
        stop("correct can only be TRUE or FALSE")
    }

    # max iteration for integration correction must be a positive integer
    if(iter_max <= 0 | !is.wholenumber(iter_max)) {
        stop("Maximum number of iterations in integration corrrection step must be a positive integer")
    }

    # check that taxrate t0 inputs are numeric
    if(!is.numeric(t0)) {
        stop("t0 must be a numeric value")
    }

    # check that taxrate t1 inputs are numeric
    if(!is.numeric(t1)) {
        stop("t1 must be a numeric value")
    }

    # flag if t0 outside unit circle
    if(t0 < 0 | t0 > 1) {
        warning("Are you sure this is the correct value for t0? Note that t0=1 means a 100% tax rate!")
    }

    # flag if t1 outside unit circle
    if(t1 < 0 | t1 > 1) {
        warning("Are you sure this is the correct value for t1? Note that t1=1 means a 100% tax rate!")
    }

    # checks of inputs t0 Vs t1
    if(t0 == t1) {
        stop("Cannot calculate elasticity (t0 cannot equal t1")
    }

    # t1 > 0
    if(t1 < t0) {
        stop("t1 must be larger than t0")
    }

    # is notch choice a logical value?
    if(!is.logical(notch)) {
        stop("notch can either be TRUE or FALSE (i.e. kink).")
    }

    # is force_notch choice a logical value?
    if(!is.logical(force_notch)) {
        stop("force_notch can either be TRUE or FALSE.")
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
        stop("p_miny must be numeric")
    }


    # is p_maxy numeric?
    if(!is.na(p_maxy) & !is.numeric(p_maxy)) {
        stop("p_maxy must be numeric (unless unspecified using NA)")
    }

    # is p_ybreaks numeric?
    if(!is.null(p_ybreaks) & !is.numeric(p_ybreaks)) {
        stop("p_ybreaks must be a numeric sequence (unless unspecified using NULL)")
    }

    # is p_theme a string?
    if(!is.character(p_theme)) {
        stop("p_theme must be a string (theme in ggplot2 format), e.g. 'bw_light")
    }

    # is p_freq_color a string?
    if(!is.character(p_freq_color)) {
        stop("p_freq_color choice must be a string, e.g. 'black'")
    }

    # is p_cf_color a string?
    if(!is.character(p_cf_color)) {
        stop("p_cf_color choice must be a string, e.g. 'maroon'")
    }

    # is p_zstar_color a string?
    if(!is.character(p_zstar_color)) {
        stop("p_zstar_color choice must be a string, e.g. 'red'")
    }


    # is p_freq_size numeric?
    if(!is.numeric(p_freq_size)) {
        stop("p_freq_size choice must be numeric")
    }

    # is p_freq_msize numeric and positive?
    if(p_freq_msize <= 0 | !is.numeric(p_freq_msize)) {
        stop("p_freq_msize must be a positive numeric value")
    }

    # is p_cf_size numeric?
    if(!is.numeric(p_cf_size)) {
        stop("p_cf_size choice must be numeric")
    }

    # is p_zstar_size numeric?
    if(!is.numeric(p_zstar_size)) {
        stop("p_zstar_size choice must be numeric")
    }

    # is p_grid_major_y_color a string?
    if(!is.character(p_grid_major_y_color)) {
        stop("p_grid_major_y_color choice must a string, e.g. 'blue'")
    }


    # is p_b a logical?
    if(!is.logical(p_b)) {
        stop("p_b (whether to show bunching estimate on plot) must be TRUE or FALSE")
    }

    # if p_b_xpos is selected by user, is p_b_xpos numeric and within the range of x?
    if(!is.na(p_b_xpos)) {
        if(!is.numeric(p_b_xpos) | p_b_xpos < data_varmin | p_b_xpos > data_varmax) {
            stop("p_b_xpos must be numeric and lie within the data's range")
        }
    }

    # if p_b_ypos is selected by user, is it numeric?
    if(!is.na(p_b_xpos)) {
        if(!is.numeric(p_b_ypos)) {
            stop("p_b_ypos must be numeric")
        }
    }

    # is p_b_size numeric and positive?
    if(p_b_size <= 0 | !is.numeric(p_b_size)) {
        stop("p_b_size must be a positive numeric value")
    }

    # is p_domregion_color a string?
    if(!is.character(p_domregion_color)) {
        stop("p_domregion_color choice must be a string, e.g. 'blue'")
    }

    # is p_domregion_ltype a string?
    if(!is.character(p_domregion_ltype)) {
        stop("p_domregion_ltype choice must be a string, e.g. 'longdash'")
    }

    # if seed not NA, is it numeric?
    if(!is.na(seed) & !is.numeric(seed)) {
        stop("seed must be numeric")
    }

    # is e_parametric logical?
    if(!is.logical(e_parametric)) {
        stop("e_parametric can only be TRUE or FALSE")
    }



    # -----------------------------
    # 1. bin the data
    # -----------------------------
    # turn data into dataframe of binned counts
    binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)

    # -----------------------------
    # 2. first pass prep and fit
    # -----------------------------

    # Kink case

    if (notch == F) {

        # prepare data
        firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                      poly, bins_excl_l, bins_excl_r, rn, extra_fe)

        # fit firstpass model
        # set zD_bin to NA if it's a kink
        zD_bin <- NA
        zD <- NA
        firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula,
                                            notch, zD_bin)

    } else if (notch == T) {
        # calculate z_dominated for notches
        z_dominated <- bunching::domregion(zstar, t0, t1, binwidth)
        zD_bin <- z_dominated$zD_bin
        zD <- z_dominated$zD

        # if we force zu, same procedure as kink
        if (force_notch == T) {
            # prepare data
            firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                          poly, bins_excl_l, bins_excl_r, rn, extra_fe)
            # fit firstpass model
            firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch, zD_bin)



        } else if (force_notch == F) {
            # start with only one bin above zstar
            bins_excl_r <- 1
            firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                          poly, bins_excl_l, bins_excl_r, rn, extra_fe)
            firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch, zD_bin)

            # extract bunching mass below and missing mass above zstar
            B_below <- firstpass$B_zl_zstar
            M_above <- -firstpass$B_zstar_zu

            # check that missing mass above is smaller. if not, stop
            if(M_above > B_below) {
                stop("Missing mass above zstar is larger than bunching mass below. Are you sure this is a notch?")
            }

            # if(M_above > B_below) == TRUE, we start shifting zu bins up until B = M.
            # first, must set upper bound on number of iterations.
            # bins above zstar cannot be more than min of:
            #   1. available bins
            #   2. remaining DoF
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
                firstpass_prep$model_formula <- stats::as.formula(paste(Reduce(paste, deparse(firstpass_prep$model_formula)), newvar, sep = " + "))
                # re-fit model using the now expanded zu
                firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch, zD_bin)
                # get new B below and M above
                B_below <- firstpass$B_zl_zstar
                M_above <- -firstpass$B_zstar_zu
            }
            # assign final zu_bin to bins_excl_r (used for plotting)
            bins_excl_r <- zu_bin
            # update bins_above_excluded to relate to this new bins_excl_r
            firstpass_prep$data_binned$bin_above_excluded <- ifelse(firstpass_prep$data_binned$z_rel > zu_bin,1,0)

        }

    }

    # after fitting is done, extract info (counterfactual, residuals, etc.)
    counterfactuals_for_graph <- firstpass$cf_density
    residuals_for_boot <- firstpass$residuals
    bunchers_initial <- firstpass$bunchers_excess
    b_estimate <- firstpass$b_estimate
    mbuncher <- bunching::marginal_buncher(beta = b_estimate, binwidth = binwidth, zstar = zstar)
    e_estimate <- bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch, e_parametric = e_parametric)
    model_fit <- firstpass$coefficients
    alpha <- firstpass$alpha




    # if we don't do correction, bunchers_initial will also be the final B_for_output
    B_for_output <- bunchers_initial

    # initialise bootstrap results
    b_sd <- NA
    b_vector <- NA
    e_vector <-NA
    e_sd <- NA
    B_vector <- NA
    B_sd <- NA
    alpha_vector <- NA
    alpha_sd <- NA
    mbuncher_vector <- NA
    mbuncher_sd <- NA

    # -----------------------------------------
    # 3. if no correction needed, do bootstrap if asked for
    # -----------------------------------------

    if(correct == F) {

        if(n_boot > 0) {
            boot_results <- bunching::do_bootstrap(zstar, binwidth, firstpass_prep, residuals_for_boot, n_boot,
                                                   correct, iter_max, notch, zD_bin, seed)
            b_sd <- boot_results$b_sd
            b_vector <- boot_results$b_vector
            e_vector <- unlist(lapply(b_vector, function(b) {
                bunching::elasticity(b, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch, e_parametric = e_parametric)
            }))
            e_sd <- round(stats::sd(e_vector),3)
            #B_for_output <- bunchers_initial # this is bunchers excess. if we dont do integration constraint, this will be output
            B_vector <- boot_results$B_vector
            B_sd <- boot_results$B_sd
            alpha_vector <- boot_results$alpha_vector
            alpha_sd <- boot_results$alpha_sd
            mbuncher_vector <- boot_results$marginal_buncher_vector
            mbuncher_sd <- boot_results$marginal_buncher_sd

        }
    }

    # -------------------------------------------------------------------------
    # 4. if correction requested, do that first to get residuals, then bootstrap
    # -------------------------------------------------------------------------
    if (correct == T) {

        # initial correction to get updated estimates, and vector of residuals for bootstrap
        firstpass_corrected <- bunching::do_correction(zstar, binwidth, firstpass_prep$data_binned,
                                                       firstpass, iter_max, notch, zD_bin)

        # get corrected results for beta, counterfactuals, alpha etc.
        counterfactuals_for_graph <- firstpass_corrected$data$cf_density
        residuals_for_boot <- firstpass_corrected$data$residuals
        B_for_output <- firstpass_corrected$B_corrected
        b_estimate <- firstpass_corrected$b_corrected
        mbuncher <- firstpass_corrected$marginal_buncher_corrected
        e_estimate <- bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch, e_parametric = e_parametric)
        model_fit <- firstpass_corrected$coefficients
        alpha <- firstpass_corrected$alpha_corrected


        # we now have the correct residuals. add to our original data in firstpass_prep$data
        # and bootstrap if requested

        if(n_boot > 0) {
            boot_results <- bunching::do_bootstrap(zstar, binwidth, firstpass_prep, residuals_for_boot, n_boot,
                                                   correct,iter_max, notch, zD_bin, seed)

            b_vector <- boot_results$b_vector
            b_sd <- boot_results$b_sd

            e_vector <- unlist(lapply(b_vector, function(b) {
                bunching::elasticity(b, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch, e_parametric = e_parametric)
            }))
            e_sd <- round(stats::sd(e_vector),3)

            B_vector <- boot_results$B_vector
            B_sd <- boot_results$B_sd
            alpha_vector <- boot_results$alpha_vector
            alpha_sd <- boot_results$alpha_sd
            mbuncher_vector <- boot_results$marginal_buncher_vector
            mbuncher_sd <- boot_results$marginal_buncher_sd
        }
    }



    # add checks for estimated values
    # b_estimate < 0 (implies elasticity < 0, and mbuncher < zstar)
    if(b_estimate < 0) {
        warning("b_estimate is less than zero. This implies negative elasticity and marginal buncher below zstar, which cannot be. \n
             Are you sure this is a notch? Check your choice of zstar, bins_excl_l and bins_excl_r")
    }

    ##################
    # notch checks
    ##################
    # is alpha within 0-1?
    if(is.numeric(alpha) & (alpha > 1 | alpha < 0)) {
        warning("The estimated alpha (fraction in dominated region) is not between 0-1. Are you sure this is a notch?")
    }


    # in case of notch, zD > mbuncher is nonsensical
    if(!is.na(zD) & (zD > mbuncher)) {
        warning("estimated zD (upper bound of dominated region) is larger than estimated marginal buncher's counterfactual z level \n Are you sure this is a notch? \n If yes, check your input choices for t0, t1, and force_notch.")
    }

    # in case of notch, zD > bins_excl_r is nonsensical
    if(!is.na(zD) & (zD > bins_excl_r)) {
        warning("estimated zD (upper bound of dominated region) is larger than bins_excl_r (upper bound of bunching region) \n Are you sure this is a notch? \n If yes, check your input choices for t0, t1, and force_notch.")
    }


    # ---------------------------------------------
    # 5. make plot
    # ---------------------------------------------

    # if p_b_xpos/p_y_xpos not chosen, set them. get data ranges
    zmin <- min(firstpass_prep$data_binned$bin)
    zmax <- max(firstpass_prep$data_binned$bin)
    maxy <- max(firstpass_prep$data_binned$freq_orig, counterfactuals_for_graph)

    # x position
    if(is.na(p_b_xpos)) {
        if (notch == T) {
            p_b_xpos <- zmin + (zstar - zmin)*.3
        } else {
            p_b_xpos <- zstar + (zmax - zstar)*.7
        }
    }

    # y position
    if(is.na(p_b_ypos)) {
        p_b_ypos <- maxy * .8
    }
    # get name of z_vector to pass as xtitle if chosen
    if (p_xtitle == "z_name") {
        p_xtitle <- deparse(substitute(z_vector))
    }

    # theme
    if (p_theme == "bw_light") {
        p_theme <- "theme_bw() + theme_light()"
    }


    # plot!
    p <- bunching::plot_bunching(firstpass_prep$data_binned, cf = counterfactuals_for_graph, zstar,
                                 binwidth, bins_excl_l, bins_excl_r,
                                 p_title, p_xtitle, p_ytitle, p_miny, p_maxy, p_ybreaks, p_axis_title_size, p_axis_val_size,
                                 p_theme, p_freq_color, p_cf_color, p_zstar_color, p_grid_major_y_color,
                                 p_freq_size, p_cf_size, p_freq_msize, p_zstar_size,
                                 p_b, b = b_estimate, b_sd = b_sd, p_e, e = e_estimate, e_sd = e_sd,
                                 p_b_xpos, p_b_ypos, p_b_size,
                                 t0, t1, notch, p_domregion_color, p_domregion_ltype, n_boot)


    output <- list("plot" = p,
                   "data" = firstpass_prep$data_binned,
                   "cf" = counterfactuals_for_graph,
                   "model_fit" = model_fit,
                   "B" = round(B_for_output,3),
                   "B_vector" = B_vector,
                   "B_sd" = B_sd,
                   "b" = b_estimate,
                   "b_vector" = b_vector,
                   "b_sd" = b_sd,
                   "e" = round(e_estimate,3),
                   "e_vector" = e_vector,
                   "e_sd" = e_sd,
                   "alpha" = round(alpha,3),
                   "alpha_vector" = alpha_vector,
                   "alpha_sd" = alpha_sd,
                   "zD" = zD,
                   "zD_bin" = zD_bin,
                   "marginal_buncher" = mbuncher,
                   "marginal_buncher_vector" = mbuncher_vector,
                   "marginal_buncher_sd" = mbuncher_sd)
    return(output)
}

