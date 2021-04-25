#' Fit Bunching
#'
#' Fit bunching model to (binned) data and estimate excess mass.
#' @param thedata (binned) data that includes all variables necessary for fitting the model.
#' @param themodelformula formula to fit.
#' @param binwidth a numeric value for the width of each bin.
#' @param zD_bin the bin marking the upper end of the dominated region (notch case).
#' @inheritParams bunchit
#' @return \code{fit_bunching} returns a list of the following results:
#' \item{coefficients}{The coefficients from the fitted model.}
#' \item{residuals}{The residuals from the fitted model.}
#' \item{cf_density}{The estimated counterfactual density.}
#' \item{bunchers_excess}{The estimate of the excess mass (not normalized).}
#' \item{cf_bunchers}{The counterfactual estimate of counts in the bunching region.}
#' \item{b_estimate}{The estimate of the normalized excess mass.}
#' \item{bins_bunchers}{The number of bins in the bunching region.}
#' \item{model_formula}{The model formula used for fitting.}
#' \item{B_zl_zstar}{The count of bunchers in the bunching region below and up to zstar.}
#' \item{B_zstar_zu}{The count of bunchers in the bunching region above zstar.}
#' \item{alpha}{The estimated fraction of bunchers in the dominated region (only in notch case.)}
#' \item{zD_bin}{The value of the bin which zD falls in.}
#' @seealso \code{\link{bunchit}}, \code{\link{prep_data_for_fit}}
#'
#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' prepped_data <- prep_data_for_fit(binned_data, zstar = 10000, binwidth = 50,
#'                                   bins_l = 20, bins_r = 20, poly = 4)
#' fitted <- fit_bunching(thedata = prepped_data$data_binned,
#'                        themodelformula = prepped_data$model_formula,
#'                        binwidth = 50)
#' # extract coefficients
#' fitted$coefficients

#' @export


fit_bunching <- function(thedata, themodelformula, binwidth, notch = FALSE, zD_bin = NA) {

    # fit model, extract coefficients and residuals (for bootstrap later)
    model_fit <- stats::lm(themodelformula,thedata)
    coefficients <- summary(model_fit)$coefficients
    residuals <- stats::residuals(model_fit)

    # estimate couterfactual
    thedata$cf <- stats::predict(model_fit, thedata)
    # remove zstar dummy effect
    thedata$cf <- thedata$cf - (thedata$zstar * coefficients["zstar", "Estimate"])

    # remove excluded region dummy effects
    bins_excluded_in_reg <- rownames(coefficients)[grepl("bin_excl",rownames(coefficients))]
    for(i in bins_excluded_in_reg){
        thedata$cf <- thedata$cf - (thedata[[i]] * coefficients[i,"Estimate"] )
    }

    # estimate bunching mass by region: outside bunching region, zl to zstar, and (above) zstar to zu
    # note: separation of zl_zstar and zstar_zu is useful for notches

    # bins_excl_r doesnt get updated here so must do it manually
    # get number of bins above zstar through formula since that will get updated for notch iterations
    bins_zstar_zu <- sum(grepl("bin_excl_r",rownames(coefficients)))
    bins_zl_zstar <- sum(grepl("bin_excl_l",rownames(coefficients))) + 1 # 1 is for the zstar point

    # get zstar value
    zstarvalue <- thedata$bin[thedata$zstar == 1]
    # use these to make indicator for zl_zstar and zstar_zu
    thedata$zl_zstar <- ifelse((thedata$bin >= zstarvalue - (binwidth * (bins_zl_zstar - 1))) & (thedata$bin <= zstarvalue), 1, 0)
    thedata$zstar_zu <- ifelse((thedata$bin <= zstarvalue + (binwidth * bins_zstar_zu)) & (thedata$bin > zstarvalue), 1, 0)

    # use these to make indicator by bunching region
    thedata$bunch_region <- ifelse(thedata$zl_zstar == 1, "zl_zstar",
                                   ifelse(thedata$zstar_zu == 1, "zstar_zu",
                                          "outside_bunching"))

    # counts by bunching_region:
    bunching_region_count <- thedata %>%
        dplyr::group_by(bunch_region) %>%
        dplyr::summarize(actual = sum(freq_orig),
                         cf = sum(cf),
                         excess = actual - cf)

    # get usual Bunching Mass estimates
    B_zl_zstar <- as.numeric(subset(bunching_region_count, bunch_region == "zl_zstar", select = "excess"))
    B_zstar_zu <- as.numeric(subset(bunching_region_count, bunch_region == "zstar_zu", select = "excess"))
    # B_zstar_zu will be NA if we have no bins excluded above zstar. set to 0
    if(is.na(B_zstar_zu)) {
        B_zstar_zu <- 0
    }
    # get B: total bunching from zl to zu (traditional estimate for kinks)
    bunchers_excess <- B_zl_zstar + B_zstar_zu

    # counterfactual bunchers
    cf_bunchers <- sum(subset(bunching_region_count, bunch_region != "outside_bunching", select = "cf"))

    # number of bins in excluded region
    bins_bunchers <- sum(thedata$bunch_region %in% c("zl_zstar", "zstar_zu"))

    # average per bin counterfactual
    c0 <- cf_bunchers/bins_bunchers

    # normalized b
    b_estimate <- as.numeric(sprintf("%.9f", bunchers_excess/c0))

    # alpha set to NA if we don't do notch below but need to pass to results
    alpha <- NA

    # --------------------------------------------------------------------
    #           if a notch, bunchers are only B_zl_zstar,
    #           bins_bunchers are only those <= zstar
    # --------------------------------------------------------------------

    if(notch) {
        bunchers_excess <- B_zl_zstar
        # number of bins in excluded region
        bins_bunchers <- sum(thedata$zl_zstar)
        # in notch case, instead of avg c0 get height of counterfactual at (zstar)
        c0 <- thedata$cf[thedata$zstar == 1]
        # normalized b
        b_estimate <- as.numeric(sprintf("%.9f", bunchers_excess/c0))
        # alpha: fraction in dominated region zstar to zD bin
        domregion_freq <- sum(thedata$freq_orig[thedata$z_rel >= 1 & thedata$z_rel <= zD_bin])
        domregion_cf <- sum(thedata$cf[thedata$z_rel >= 1 & thedata$z_rel <= zD_bin])
        alpha <- domregion_freq/domregion_cf

    }

    # return output we need
    output <- list("coefficients" = coefficients,
                   "residuals" = residuals,
                   "cf_density" = thedata$cf,
                   "c0" = c0,
                   "bunchers_excess" = bunchers_excess,
                   "cf_bunchers" = cf_bunchers,
                   "b_estimate" = b_estimate,
                   "bins_bunchers" = bins_bunchers,
                   "model_formula" = themodelformula,
                   "B_zl_zstar" = B_zl_zstar,
                   "B_zstar_zu" = B_zstar_zu,
                   "alpha" = alpha,
                   "zD_bin" = zD_bin)

    return(output)
}
