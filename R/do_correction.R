#' Integration Constraint Correction
#'
#' Implements the correction for the integration constraint.
#'
#' @param data_prepped (binned) data that includes all variables necessary for fitting the model.
#' @param firstpass_results initial bunching estimates without correction.
#' @inheritParams bunchit
#' @inheritParams fit_bunching
#' @seealso \code{\link{bunchit}}, \code{\link{fit_bunching}}
#' @return do_correction returns a list with the data and estimates after correcting for the integration constraint, as follows:
#' \item{data}{The dataset with the corrected counterfactual.}
#' \item{coefficients}{The coefficients of the model fit on the corrected data.}
#' \item{b_corrected}{The normalized excess mass, corrected for the integration constraint.}
#' \item{B_corrected}{The excess mass (not normalized), corrected for the integration constraint.}
#' \item{c0_corrected}{The counterfactual at zstar, corrected for the integration constraint.}
#' \item{marginal_buncher_corrected}{The location (z value) of the marginal buncher, corrected for the integration constraint.}
#' \item{alpha_corrected}{The estimated fraction of bunchers in the dominated region, corrected for the integration constraint (only in notch case).}

#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' prepped_data <- prep_data_for_fit(binned_data, zstar = 10000, binwidth = 50,
#'                                   bins_l = 20, bins_r = 20, poly = 4)
#' firstpass <- fit_bunching(prepped_data$data_binned, prepped_data$model_formula)
#' corrected <- do_correction(zstar = 10000, binwidth = 50,
#'                            data_prepped = prepped_data$data_binned,
#'                            firstpass_results = firstpass)
#' paste0("Without correction, b = ", firstpass$b_estimate)
#' paste0("With correction, b = ", round(corrected$b_corrected,3))
#' @export



do_correction <- function(zstar, binwidth, data_prepped, firstpass_results,
                          correct_iter_max = 200, notch = FALSE, zD_bin = NA) {
    # get initial buncher value, bins of bunchers and model formula
    bunchers_excess_initial <- firstpass_results$bunchers_excess
    bins_bunchers <- firstpass_results$bins_bunchers
    model_formula <- firstpass_results$model_formula

    # calculate proportional shift upwards for those above zU
    data_prepped$location_shift_sca <- data_prepped$bin_above_excluded/
        sum(data_prepped$freq[data_prepped$bin_above_excluded == 1]) # total count above excluded region


    # set some initial value
    b_diff <- 1000

    # create df for estimates of iterations
    excess_updated_df <- data.frame("iteration" = 0,"bunchers_excess_updated" = bunchers_excess_initial)
    bunchers_excess_updated <- bunchers_excess_initial # to pass as variable below

    iteration <- 1
    while(b_diff >= 1 & iteration < correct_iter_max){
        data_prepped$freq <- data_prepped$freq_orig * (1 + (bunchers_excess_updated*data_prepped$location_shift_sca))
        iteration_results <- bunching::fit_bunching(data_prepped, model_formula, notch, zD_bin)
        bunchers_excess_updated <- iteration_results$bunchers_excess
        c0_updated <- iteration_results$c0
        # add data to excess_updated_df
        res <- c(as.integer(iteration),bunchers_excess_updated)
        excess_updated_df <- rbind(excess_updated_df,res)

        #  first iteration is 0 (the original count)
        b_diff <- excess_updated_df$bunchers_excess_updated[iteration] -
            excess_updated_df$bunchers_excess_updated[iteration+1]
        iteration <- iteration + 1
    }

    # b_corrected
    b_corrected <- bunchers_excess_updated/c0_updated

    # B corrected (excess mass without normalization)
    B_corrected <- bunchers_excess_updated
    # if we did correction, get last cf (this is the correct one) and assign to cf_graph for graphing
    data_prepped$cf_density <- iteration_results$cf_density

    # get alpha
    alpha_corrected <- iteration_results$alpha

    # get marginal buncher
    mbuncher_corrected <- bunching::marginal_buncher(beta = b_corrected, binwidth = binwidth, zstar = zstar,
                                                     notch = notch, alpha = alpha_corrected)
    # get residuals
    data_prepped$residuals <- data_prepped$cf_density - data_prepped$freq_orig



    # generate output (updated with correction)
    output <- list("data" = data_prepped,
                   "coefficients" = iteration_results$coefficients,
                   "b_corrected" = b_corrected,
                   "B_corrected" = B_corrected,
                   "c0_corrected" = c0_updated,
                   "marginal_buncher_corrected" = mbuncher_corrected,
                   "alpha_corrected" = alpha_corrected)

    return(output)
}
