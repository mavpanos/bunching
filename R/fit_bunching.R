#' Fit  bunching model and estimate excess mass.
#' @param thedata (binned) data that includes all variables necessary for fitting the model.
#' @param themodelformula formula to fit.
#' @param zD_bin the bin marking the upper end of the dominated region (notch case).
#' @inheritParams bunchit
#' @return coefficients, residuals, cf_density, bunchers_excess, cf_bunchers, b_estimate, bins_bunchers, model_formula.
#' @seealso \code{\link{bunchit}}, \code{\link{prep_data_for_fit}}
#' @export


fit_bunching <- function(thedata, themodelformula, notch, zD_bin) {

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
    # separation of zl_zstar and zstar_zu is useful for notches

    # bins_excl_r doesnt get updated here so must do it manually
    # get number of bins above zstar through formula since that will get updated for notch iterations
    bins_zstar_zu <- sum(grepl("bin_excl_r",rownames(coefficients)))
    bins_zl_zstar <- sum(grepl("bin_excl_l",rownames(coefficients))) + 1 # 1 is for the zstar point

    # get zstar value
    zstarvalue <- thedata$bin[thedata$zstar == 1]
    # get binwidth value (difference between binvalue of consecutive rows)
    binwidthvalue <- thedata$bin[2] - thedata$bin[1]
    # use these to make indicator for zl_zstar and zstar_zu
    thedata$zl_zstar <- ifelse((thedata$bin >= zstarvalue - (binwidthvalue * (bins_zl_zstar - 1))) & (thedata$bin <= zstarvalue), 1, 0)
    thedata$zstar_zu <- ifelse((thedata$bin <= zstarvalue + (binwidthvalue * bins_zstar_zu)) & (thedata$bin > zstarvalue), 1, 0)

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

    # normalised b
    b_estimate <- as.numeric(sprintf("%.3f", bunchers_excess/c0))

    # alpha set to NA if we don't do notch below but need to pass to results
    alpha <- NA
    ################################################################################
    # if Notch, bunchers are only B_zl_zstar, bins_bunchers are only those <= zstar
    ################################################################################
    if(notch) {
        bunchers_excess <- B_zl_zstar
        # number of bins in excluded region
        bins_bunchers <- sum(thedata$zl_zstar)
        # average per bin counterfactual
        c0 <- cf_bunchers/bins_bunchers
        # normalised b
        b_estimate <- as.numeric(sprintf("%.3f", bunchers_excess/c0))
        # alpha: fraction in dominated region zstar to zD bin
        domregion_freq <- sum(thedata$freq_orig[thedata$z_rel >= 1 & thedata$z_rel <= zD_bin])
        domregion_cf <- sum(thedata$cf[thedata$z_rel >= 1 & thedata$z_rel <= zD_bin])
        alpha <- domregion_freq/domregion_cf

    }



    # return output we need
    output <- list("coefficients" = coefficients,
                   "residuals" = residuals,
                   "cf_density" = thedata$cf,
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
