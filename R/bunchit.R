#' Bunching Estimator
#'
#' Implement the bunching estimator in a kink or notch setting.
#'
#' @param z_vector a numeric vector of (unbinned) data.
#' @param binv a string setting location of zstar within its bin ("min", "max" or "median" value). Default is median.
#' @param zstar a numeric value for the the bunching point.
#' @param binwidth a numeric value for the width of each bin.
#' @param bins_l number of bins to left of zstar to use in analysis.
#' @param bins_r number of bins to right of zstar to use in analysis.
#' @param poly a numeric value for the order of polynomial for counterfactual fit. Default is 9.
#' @param bins_excl_l number of bins to left of zstar to include in bunching region. Default is 0.
#' @param bins_excl_r number of bins to right of zstar to include in bunching region. Default is 0.
#' @param extra_fe a numeric vector of bin values to control for using fixed effects. Default includes no controls.
#' @param rn a numeric vector of (up to 2) round numbers to control for. Default includes no controls.
#' @param n_boot number of bootstrapped iterations. Default is 100.
#' @param correct implements correction for integration constraint. Default is TRUE.
#' @param correct_above_zu if integration constraint correction is implemented, should counterfactual be shifted only above zu (upper bound of exclusion region)? Default is FALSE (i.e. shift from above zstar).
#' @param correct_iter_max maximum iterations for integration constraint correction. Default is 200.
#' @param seed a numeric value for bootstrap seed (random re-sampling of residuals). Default is NA.
#' @param p_title plot's title. Default is empty.
#' @param p_xtitle plot's x_axis label. Default is the name of z_vector.
#' @param p_ytitle plot's y_axis label. Default is "Count".
#' @param p_miny plot's minimum y_axis value. Default is 0.
#' @param p_maxy plot's maximum y_axis value. Default is optimized internally.
#' @param p_ybreaks a numeric vector of y-axis values at which to add horizontal line markers in plot. Default is optimized internally.
#' @param p_title_size size of plot's title. Default is 11.
#' @param p_axis_title_size size of plot's axes' title labels. Default is 10.
#' @param p_axis_val_size size of plot's axes' numeric labels. Default is 8.5.
#' @param p_freq_color plot's frequency line color. Default is "black".
#' @param p_cf_color plot's counterfactual line color. Default is "maroon".
#' @param p_zstar_color plot's bunching region marker lines color. Default is "red".
#' @param p_grid_major_y_color plot's y-axis major grid line color. Default is "lightgrey".
#' @param p_freq_size plot's frequency line thickness. Default is 0.5.
#' @param p_cf_size plot's counterfactual line thickness. Default is 0.5.
#' @param p_freq_msize plot's frequency line marker size. Default is 1.
#' @param p_zstar_size plot's bunching region marker line thickness. Default is 0.5.
#' @param p_b whether plot should also include the bunching estimate. Default is FALSE.
#' @param p_e whether plot should also include the elasticity estimate. Only shown if p_b is TRUE. Default is FALSE.
#' @param p_b_e_xpos plot's x-axis coordinate of bunching/elasticity estimate. Default is set internally.
#' @param p_b_e_ypos plot's y-axis coordinate of bunching/elasticity estimate. Default is set internally.
#' @param p_b_e_size size of plot's printed bunching/elasticity estimate. Default is 3.
#' @param t0 numeric value setting the marginal (average) tax rate below zstar in a kink (notch) setting.
#' @param t1 numeric value setting the marginal (average) tax rate above zstar in a kink (notch) setting.
#' @param notch whether analysis is for a kink or notch. Default is FALSE (kink).
#' @param force_notch whether to enforce user's choice of zu (upper limit of bunching region) in a notch setting. Default is FALSE (zu set by setting bunching equal to missing mass).
#' @param e_parametric whether to estimate elasticity using parametric specification (quasi-linear and iso-elastic utility function). Default is FALSE (which estimates reduced-form approximation).
#' @param e_parametric_lb lower bound for elasticity estimate's solution using parametric specification in notch setting. Default is 1e-04.
#' @param e_parametric_ub upper bound for elasticity estimate's solution using parametric specification in notch setting. Default is 3.
#' @param p_domregion_color plot's dominated region marker line color in notch setting. Default is "blue".
#' @param p_domregion_ltype line type for the vertical line type marking the dominated region (zD) in the plot for notch settings. Default is "longdash".

#' @details bunchit implements the bunching estimator in both kink and notch settings. It bins a given numeric vector, fits a counterfactual density, and estimates the bunching mass (normalized and not), the elasticity and the location of the marginal buncher. In the case of notches, it also finds the dominated region and estimates the fraction of observations located in it.

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
#'   \item{zU_bin}{The location of zU (upper range of excluded region) as estimated from notch setting by setting force_notch = FALSE.}
#'   \item{marginal_buncher}{The location (z value) of the marginal buncher.}
#'   \item{marginal_buncher_vector}{The vector of bootstrapped marginal_buncher values.}
#'   \item{marginal_buncher_sd}{The standard deviation of marginal_buncher_vector.}
#' @import BB
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @seealso \code{\link{plot_hist}}
#'
#' @examples
#' # First, load the example data
#' data(bunching_data)
#'
#' # Example 1: Kink with integration constraint correction
#' kink1 <- bunchit(z_vector = bunching_data$kink, zstar = 10000, binwidth = 50,
#'                  bins_l = 20, bins_r = 20, poly = 4, t0 = 0, t1 = .2,
#'                  p_b = TRUE, seed = 1)
#' kink1$plot
#' kink1$b
#' kink1$b_sd
#'
#' # Example 2: Kink with diffuse bunching
#' bpoint <- 10000; binwidth <- 50
#' kink2_vector <- c(bunching_data$kink_vector,
#'                  rep(bpoint - binwidth,80), rep(bpoint - 2*binwidth,190),
#'                  rep(bpoint + binwidth,80), rep(bpoint + 2*binwidth,80))
#' kink2 <- bunchit(z_vector = kink2_vector, zstar = 10000, binwidth = 50,
#'                  bins_l = 20, bins_r = 20, poly = 4,  t0 = 0, t1 = .2,
#'                  bins_excl_l = 2, bins_excl_r = 2, correct = FALSE,
#'                  p_b = TRUE, seed = 1)
#' kink2$plot
#'
#' # Example 3: Kink with further bunching at other level in bandwidth
#' kink3_vector <- c(bunching_data$kink_vector, rep(10200,540))
#' kink3 <- bunchit(kink3_vector, zstar = 10000, binwidth = 50,
#'                  bins_l = 40, bins_r = 40, poly = 6, t0 = 0, t1 = .2,
#'                  correct = FALSE, p_b = TRUE, extra_fe = 10200, seed = 1)
#' kink3$plot
#'
#' # Example 4: Kink with round number bunching
#' rn1 <- 500;  rn2 <- 250
#' bpoint <- 10000
#' kink4_vector <- c(bunching_data$kink_vector,
#'                   rep(bpoint + rn1, 270),
#'                   rep(bpoint + 2*rn1,230),
#'                   rep(bpoint - rn1,260),
#'                   rep(bpoint - 2*rn1,275),
#'                   rep(bpoint + rn2, 130),
#'                   rep(bpoint + 3*rn2,140),
#'                   rep(bpoint - rn2,120),
#'                   rep(bpoint - 3*rn2,135))
#' kink4 <- bunchit(z_vector = kink4_vector, zstar = bpoint, binwidth = 50,
#'                  bins_l = 20, bins_r = 20, poly = 6, t0 = 0, t1 = .2,
#'                  correct = FALSE, p_b = TRUE, p_e = TRUE, p_freq_msize = 1.5,
#'                  p_b_e_ypos = 880, rn = c(250,500), seed = 1)
#' kink4$plot
#'
#' # Example 5: Notch
#' notch <- bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
#'                  bins_l = 40, bins_r = 40, poly = 5, t0 = 0.18, t1 = .25,
#'                  correct = FALSE, notch = TRUE,p_b = TRUE, p_b_e_xpos = 8900,
#'                  n_boot = 0)
#' notch$plot


#' @export

bunchit <- function(z_vector, binv = "median", zstar, binwidth, bins_l, bins_r,
                    poly = 9, bins_excl_l = 0, bins_excl_r = 0, extra_fe = NA, rn = NA,
                    n_boot = 100, correct = TRUE, correct_above_zu = FALSE, correct_iter_max = 200,
                    t0, t1, notch = FALSE, force_notch = FALSE, e_parametric = FALSE,
                    e_parametric_lb = 0.0001, e_parametric_ub = 3, seed = NA,
                    p_title = "", p_xtitle = deparse(substitute(z_vector)), p_ytitle = "Count", p_title_size = 11,
                    p_axis_title_size = 10, p_axis_val_size = 8.5, p_miny = 0, p_maxy = NA, p_ybreaks = NA,
                    p_freq_color = "black", p_cf_color = "maroon", p_zstar_color = "red", p_grid_major_y_color = "lightgrey",
                    p_freq_size = .5, p_freq_msize = 1, p_cf_size = .5, p_zstar_size = .5,
                    p_b = FALSE, p_e = FALSE, p_b_e_xpos = NA, p_b_e_ypos = NA, p_b_e_size = 3,
                    p_domregion_color = "blue", p_domregion_ltype="longdash") {

    # ------------------------------------------------
    #           0. check inputs
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

    # zstar cannot be zero (elasticity needs Dz/zstar estimate)
    if(zstar == 0) {
        stop("zstar cannot be zero. If this your true bunching point, must re-centre it away from zero")
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

    # polynomial order cannot be negative or not a whole number
    if(poly < 0 | !is.wholenumber(poly)) {
        stop("poly must be a non-negative integer")
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
        stop("Bunching region below zstar too wide. Increase bins_l relative to bins_excl_l")
    }

    # is bunching region above zstar too wide?
    if(bins_excl_r >= bins_r - 5) {
        stop("Bunching region above zstar too wide. Increase bins_r relative to bins_excl_r")
    }

    # check that rn does not include 0
    if( 0 %in% rn) {
        stop("rn cannot include zero as a round number")
    }

    # if not NA, are all round numbers integers?
    if(sum(is.na(rn)) == 0 & sum(is.wholenumber(rn)) != length(rn)) {
        stop("Round number(s) must be integer(s)")
    }

    # check that no more than two round numbers are given
    if(length(rn) > 2) {
        stop("rn cannot include more than two unique levels for round number bunching")
    }

    # if two round numbers given, check that they are unique
    if(length(rn) == 2 & length(unique(rn)) != 2) {
        stop("the two round numbers in rn cannot be identical")
    }

    # if rn is not NA, are round numbers included in data range?
    if(sum(is.na(rn)) == 0 & sum(rn > data_varmax) != 0) {
        stop("rn includes round numbers outside of z_vector's range of values")
    }

    # if rn is not NA,
    # is data range too small to allow for round numbers? (e.g. data range is 0-400 but rn = 500)
    if(sum(is.na(rn)) == 0 & sum(rn > (data_varmax-data_varmin))) {
        stop("rn includes round numbers that are too large for z_vector's range of values. \n Use smaller values or increase bins_l and bins_r")
    }

    # number of bootstrap iterations must be non-negative integer
    if(n_boot < 0 | !is.wholenumber(n_boot)) {
        stop("bootstrap sample size n_boot must be a non-negative integer")
    }

    # flag if bootstrap samples less than 100
    if(n_boot > 0 & n_boot < 100) {
        warning(paste0("Are you sure your choice for the bootstrap sample size n_boot is large enough?"))
    }


    # correct must be logical
    if(!is.logical(correct)) {
        stop("correct (for integration constraint) can only be TRUE or FALSE")
    }

    # correct_above_zu must be logical
    if(!is.logical(correct_above_zu)) {
        stop("correct_above_zu (for integration constraint in notch setting) can only be TRUE or FALSE")
    }

    # max iteration for integration correction must be a positive integer
    if(correct_iter_max <= 0 | !is.wholenumber(correct_iter_max)) {
        stop("Maximum number of iterations (for integration constraint) correct_iter_max must be a positive integer")
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
        warning("Are you sure this is the correct value for t0? \n Note that t0 = 0 means a 0% and t0 = 1 a 100% tax rate!")
    }

    # flag if t1 outside unit circle
    if(t1 < 0 | t1 > 1) {
        warning("Are you sure this is the correct value for t1? \n Note that t1 = 0 means a 0% and t1 = 1 a 100% tax rate!")
    }

    # are t0 and t1 different?
    if(t0 == t1) {
        stop("Cannot calculate elasticity (t0 cannot equal t1")
    }

    # t1 > 0
    if(t1 < t0) {
        stop("t1 must be larger than t0")
    }

    # is notch choice a logical value?
    if(!is.logical(notch)) {
        stop("notch can either be TRUE or FALSE (i.e. kink)")
    }

    # is force_notch choice a logical value?
    if(!is.logical(force_notch)) {
        stop("force_notch can either be TRUE or FALSE")
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

    # is p_freq_color a string?
    if(!is.character(p_freq_color)) {
        stop("p_freq_color must be a string, e.g. 'black'")
    }

    # is p_cf_color a string?
    if(!is.character(p_cf_color)) {
        stop("p_cf_color must be a string, e.g. 'maroon'")
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

    # is p_cf_size numeric?
    if(p_cf_size <= 0 | !is.numeric(p_cf_size)) {
        stop("p_cf_size must be a positive numeric value")
    }

    # is p_zstar_size numeric?
    if(p_zstar_size <= 0 | !is.numeric(p_zstar_size)) {
        stop("p_zstar_size must be a positive numeric value")
    }

    # is p_grid_major_y_color a string?
    if(!is.character(p_grid_major_y_color)) {
        stop("p_grid_major_y_color must a string, e.g. 'blue'")
    }

    # is p_b a logical?
    if(!is.logical(p_b)) {
        stop("p_b (whether to show bunching estimate on plot) must be TRUE or FALSE")
    }

    # if p_b_e_xpos is selected by user, is p_b_e_xpos numeric and within the range of x?
    if(!is.na(p_b_e_xpos)) {
        if(!is.numeric(p_b_e_xpos) | p_b_e_xpos < data_varmin | p_b_e_xpos > data_varmax) {
            stop("p_b_e_xpos must be numeric and lie within z_vector's range")
        }
    }

    # if p_b_e_ypos is selected by user, is it numeric?
    if(!is.na(p_b_e_ypos)) {
        if(!is.numeric(p_b_e_ypos)) {
            stop("p_b_e_ypos must be numeric")
        }
    }

    # is p_b_e_size numeric and positive?
    if(p_b_e_size <= 0 | !is.numeric(p_b_e_size)) {
        stop("p_b_e_size must be a positive numeric value")
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

    # is e_parametric_lb numeric?
    if(!is.numeric(e_parametric_lb)) {
        stop("e_parametric_lb must be numeric")
    }

    # is e_parametric_ub numeric?
    if(!is.numeric(e_parametric_ub)) {
        stop("e_parametric_ub must be numeric")
    }

    # is e_parametric_ub > e_parametric_lb?
    if(e_parametric_ub <= e_parametric_lb) {
        stop("e_parametric_ub must be larger than e_parametric_lb")
    }


    # ------------------------------------------------
    #               1. bin the data
    # ------------------------------------------------
    # convert vector into dataframe of binned counts
    binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)

    # --------------------------------------------
    # checks: if extra kinks do not correspond
    #           to a bin, flag it
    # --------------------------------------------
    if(sum(!is.na(extra_fe) > 0) & (length(extra_fe) != sum(extra_fe %in% binned_data$bin))) {
        warning("Warning: extra_fe do not correspond to any bin value. \n Consider using the plot_hist() function and check the returned data's bin values")
    }

    # ------------------------------------------------
    #           2. first pass prep and fit
    # ------------------------------------------------

    # Kink case

    if (notch == FALSE) {

        # prepare data
        firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                      poly, bins_excl_l, bins_excl_r, rn, extra_fe,
                                                      correct_above_zu)

        # fit firstpass model
        # set alpha, zD_bin to NA if it's a kink
        alpha <- NA
        zD_bin <- NA
        zD <- NA
        zU_notch <- NA
        firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula,
                                            notch, zD_bin)

    } else if (notch == TRUE) {
        # calculate z_dominated for notches
        z_dominated <- bunching::domregion(zstar, t0, t1, binwidth)
        zD_bin <- z_dominated$zD_bin
        zD <- z_dominated$zD
        zU_notch <- bins_excl_r

        # if we force zu, same procedure as kink
        if (force_notch == TRUE) {
            # prepare data
            firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                          poly, bins_excl_l, bins_excl_r, rn, extra_fe,
                                                          correct_above_zu)
            # fit firstpass model
            firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula, notch, zD_bin)



        } else if (force_notch == FALSE) {
            # start with only one bin above zstar
            bins_excl_r <- 1
            firstpass_prep <- bunching::prep_data_for_fit(binned_data, zstar, binwidth, bins_l, bins_r,
                                                          poly, bins_excl_l, bins_excl_r, rn, extra_fe,
                                                          correct_above_zu)
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

            # create temporary object for firstpass_prep (we will be editing this as we expand the window for zu_bin)
            tmp_firstpass_prep <- firstpass_prep

            while ((B_below > M_above) & (zu_bin < notch_iterations_bound)) {
                # keep previous zu_bin in case next gives us B > M
                zu_bin_final <- zu_bin
                # now try expanding by 1 bin
                zu_bin <- zu_bin + 1
                # add next bin_excl_r dummy as column to data
                newvar <- paste0("bin_excl_r_", zu_bin)
                tmp_firstpass_prep$data_binned[[newvar]]  <- ifelse(tmp_firstpass_prep$data_binned$z_rel == zu_bin,1,0)
                # add next order bin_excl_r to formula
                tmp_firstpass_prep$model_formula <- stats::as.formula(paste(Reduce(paste, deparse(tmp_firstpass_prep$model_formula)), newvar, sep = " + "))
                # re-fit model using the now expanded zu
                tmp_firstpass <- bunching::fit_bunching(tmp_firstpass_prep$data_binned, tmp_firstpass_prep$model_formula, notch, zD_bin)
                # get new B below and M above
                B_below <- tmp_firstpass$B_zl_zstar
                M_above <- -tmp_firstpass$B_zstar_zu

                # if M is not larger, update firstpass with this (last iterations of loop will have M > B)
                if(B_below > M_above) {
                    firstpass_prep <- tmp_firstpass_prep
                    firstpass <- tmp_firstpass
                    zu_bin_final <- zu_bin
                }

            }
            # assign final zu_bin_final to bins_excl_r (used for plotting)
            bins_excl_r <- zu_bin_final
            # update zU_notch with this (to output)
            zU_notch <- bins_excl_r
            # if we chose correct_above_zu == T, i.e. shift only above zu,
            # update bins_above_excluded to relate to this new bins_excl_r = zu_bin_final
            if(correct_above_zu) {
                firstpass_prep$data_binned$bin_above_excluded <- ifelse(firstpass_prep$data_binned$z_rel > zu_bin_final,1,0)
            }
        }

    }

    # after fitting is done, extract info (counterfactual, residuals, etc.)
    counterfactuals_for_graph <- firstpass$cf_density
    residuals_for_boot <- firstpass$residuals
    bunchers_initial <- firstpass$bunchers_excess
    b_estimate <- firstpass$b_estimate
    e_estimate <- bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch,
                                       e_parametric = e_parametric, e_parametric_lb = e_parametric_lb, e_parametric_ub = e_parametric_ub)
    model_fit <- firstpass$coefficients
    alpha <- firstpass$alpha
    mbuncher <- bunching::marginal_buncher(beta = b_estimate, binwidth = binwidth, zstar = zstar, notch = notch, alpha = alpha)

    # if we don't do correction, bunchers_initial will also be the final B_for_output
    B_for_output <- bunchers_initial

    # initialize bootstrap results
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

    # ------------------------------------------------
    #   3. if no correction needed, do bootstrap
    #               if requested
    # ------------------------------------------------

    if(correct == FALSE) {

        if(n_boot > 0) {
            boot_results <- bunching::do_bootstrap(zstar, binwidth, firstpass_prep, residuals_for_boot, n_boot,
                                                   correct, correct_iter_max, notch, zD_bin, seed)
            b_sd <- boot_results$b_sd
            b_vector <- boot_results$b_vector
            e_vector <- unlist(lapply(b_vector, function(b) {
                bunching::elasticity(b, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch,
                                     e_parametric = e_parametric, e_parametric_lb = e_parametric_lb, e_parametric_ub = e_parametric_ub)
            }))
            e_sd <- stats::sd(e_vector)
            B_vector <- boot_results$B_vector
            B_sd <- boot_results$B_sd
            alpha_vector <- boot_results$alpha_vector
            alpha_sd <- boot_results$alpha_sd
            mbuncher_vector <- boot_results$marginal_buncher_vector
            mbuncher_sd <- boot_results$marginal_buncher_sd

        }
    }

    # ------------------------------------------------
    # 4. if correction requested, do that first
    #       to get residuals, then bootstrap
    # ------------------------------------------------
    if (correct == TRUE) {

        # initial correction to get updated estimates, and vector of residuals for bootstrap
        firstpass_corrected <- bunching::do_correction(zstar, binwidth, firstpass_prep$data_binned,
                                                       firstpass, correct_iter_max, notch, zD_bin)

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
                                                   correct, correct_iter_max, notch, zD_bin, seed)

            b_vector <- boot_results$b_vector
            b_sd <- boot_results$b_sd

            e_vector <- unlist(lapply(b_vector, function(b) {
                bunching::elasticity(b, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1, notch = notch,
                                     e_parametric = e_parametric,  e_parametric_lb = e_parametric_ub, e_parametric_ub = e_parametric_ub)
            }))
            e_sd <- stats::sd(e_vector)

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

    # ------------------------------------------------
    #   5. do some checks for notches' estimates
    # ------------------------------------------------

    # is alpha within 0-1?
    if(is.numeric(alpha)) {
        if(alpha > 1 | alpha < 0) {
            warning("The estimated alpha (fraction in dominated region) is not between 0-1. Are you sure this is a notch?")
        }
    }


    # in case of notch, zD > mbuncher is nonsensical
    if(!is.na(zD) & (zD > mbuncher)) {
        warning("estimated zD (upper bound of dominated region) is larger than estimated marginal buncher's counterfactual z level \n Are you sure this is a notch? \n If yes, check your input choices for t0, t1, force_notch, correct and correct_above_zu.")
    }

    # in case of notch, zD_bin > bins_excl_r is nonsensical
    if(!is.na(zD_bin) & (zD_bin > bins_excl_r)) {
        warning("estimated zD (upper bound of dominated region) is larger than bins_excl_r (upper bound of bunching region) \n Are you sure this is a notch? \n If yes, check your input choices for t0, t1, force_notch, correct and correct_above_zu.")
    }


    # ---------------------------------------------
    #               6. make plot
    # ---------------------------------------------


    # plot!
    #
    p <- bunching::plot_bunching(z_vector, binned_data = firstpass_prep$data_binned, cf = counterfactuals_for_graph, zstar = zstar,
                                 binwidth = binwidth, bins_excl_l = bins_excl_l, bins_excl_r = bins_excl_r,
                                 p_title = p_title, p_xtitle = p_xtitle, p_ytitle = p_ytitle, p_miny = p_miny,
                                 p_maxy = p_maxy, p_ybreaks = p_ybreaks, p_title_size = p_title_size,
                                 p_axis_title_size = p_axis_title_size, p_axis_val_size = p_axis_val_size,
                                 p_freq_color = p_freq_color, p_cf_color = p_cf_color, p_zstar_color = p_zstar_color,
                                 p_grid_major_y_color = p_grid_major_y_color, p_freq_size = p_freq_size,
                                 p_freq_msize = p_freq_msize, p_cf_size = p_cf_size, p_zstar_size = p_zstar_size,
                                 p_b = p_b, b = b_estimate, b_sd = b_sd, p_e = p_e, e = e_estimate, e_sd = e_sd,
                                 p_b_e_xpos = p_b_e_xpos, p_b_e_ypos = p_b_e_ypos, p_b_e_size = p_b_e_size,
                                 t0 = t0, t1 = t1, notch = notch, p_domregion_color = p_domregion_color,
                                 p_domregion_ltype = p_domregion_ltype)
    # set rounding
    round_dp <- 3

    output <- list("plot" = p,
                   "data" = firstpass_prep$data_binned,
                   "cf" = counterfactuals_for_graph,
                   "model_fit" = model_fit,
                   "B" = round(B_for_output,round_dp),
                   "B_vector" = round(B_vector,round_dp),
                   "B_sd" = round(B_sd,round_dp),
                   "b" = round(b_estimate,round_dp),
                   "b_vector" = round(b_vector,round_dp),
                   "b_sd" = round(b_sd,round_dp),
                   "e" = round(e_estimate,round_dp),
                   "e_vector" = round(e_vector,round_dp),
                   "e_sd" = round(e_sd,round_dp),
                   "alpha" = round(alpha,round_dp),
                   "alpha_vector" = alpha_vector,
                   "alpha_sd" = round(alpha_sd,round_dp),
                   "zD" = round(zD,round_dp),
                   "zD_bin" = zD_bin,
                   "zU_bin" = zU_notch,
                   "marginal_buncher" = round(mbuncher,round_dp),
                   "marginal_buncher_vector" = round(mbuncher_vector,round_dp),
                   "marginal_buncher_sd" = round(mbuncher_sd,round_dp))
    return(output)
}

