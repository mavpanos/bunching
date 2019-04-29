#' Fit  bunching model and estimate excess mass.
#' @param thedata (binned) data that includes all variables necessary for fitting the model.
#' @param themodelformula formula to fit.
#' @return coefficients, residuals, cf_density, bunchers_excess, cf_bunchers, b_estimate, bins_bunchers, model_formula.
#' @seealso \code{\link{bunchit}}, \code{\link{prepare_data}}
#' @export


fit_bunching <- function(thedata, themodelformula) {
    model_fit <- stats::lm(themodelformula,thedata)
    coefficients <- summary(model_fit)$coefficients
    residuals <- stats::residuals(model_fit)

    # estimate couterfactual
    thedata$cf <- stats::predict(model_fit, thedata)
    # remove kink dummy effect
    thedata$cf <- thedata$cf - (thedata$kink * coefficients["kink", "Estimate"])

    # remove excluded region dummy effects
    bins_excluded_in_reg <- rownames(coefficients)[grepl("bin_excl",rownames(coefficients))]
    for(i in bins_excluded_in_reg){
        thedata$cf <- thedata$cf - (thedata[[i]] * coefficients[i,"Estimate"] )
    }

    # number of bins in excluded region
    bins_bunchers <- sum(thedata$bunch_region)

    # counts by bunching_region:
    bunching_region_count <- thedata %>%
        dplyr::group_by(bunch_region) %>%
        dplyr::summarize(actual = sum(freq_orig),
                         cf = sum(cf),
                         excess = actual - cf)

    # find c0 (avg per bin counterfactual count in excluded region)
    bunchers_excess <- as.numeric(bunching_region_count[which(bunching_region_count$bunch_region == 1),"excess"])
    cf_bunchers <- as.numeric(bunching_region_count[which(bunching_region_count$bunch_region == 1),"cf"])
    c0 <- cf_bunchers/bins_bunchers

    # get b
    b_estimate <- as.numeric(sprintf("%.3f", bunchers_excess/c0))

    # return output we need
    output <- list("coefficients" = coefficients,
                   "residuals" = residuals,
                   "cf_density" = thedata$cf,
                   "bunchers_excess" = bunchers_excess,
                   "cf_bunchers" = cf_bunchers,
                   "b_estimate" = b_estimate,
                   "bins_bunchers" = bins_bunchers,
                   "model_formula" = themodelformula)
    return(output)
}
