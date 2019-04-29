#' Implements correction for integration constraint.
#'
#' @param thedata (binned) data that includes all variables necessary for fitting the model.
#' @param firstpass_results initial bunching estimates without correction.
#' @param max_iterations maximum number of iterations for counterfactual shifting.
#' @seealso \code{\link{bunchit}} \code{\link{fit_bunching}}
#' @return The data used for estimation, and the "corrected" estimates of both the pure excess mass (B_corrected) and the normalised excess mass (b_corrected).
#' @export



do_correction <- function(thedata, firstpass_results, max_iterations) {

    bunchers_excess_initial <- firstpass_results$bunchers_excess
    bins_bunchers <- firstpass_results$bins_bunchers
    model_formula <- firstpass_results$model_formula

    thedata$location_shift_sca <- thedata$bin_above_excluded/
        sum(thedata$freq[thedata$bin_above_excluded == 1]) # total count above excluded region


    # set some initial value
    b_diff <- 1000

    # create df for estimates of iterations
    excess_updated_df <- data.frame("iteration" = 0,"bunchers_excess_updated" = bunchers_excess_initial)
    bunchers_excess_updated <- bunchers_excess_initial # to pass as variable below

    iteration <- 1
    while(b_diff >= 1 & iteration < max_iterations){
        thedata$freq <- thedata$freq_orig * (1 + (bunchers_excess_updated*thedata$location_shift_sca))
        iteration_results <- run_reg(thedata, model_formula)
        bunchers_excess_updated <- iteration_results$bunchers_excess
        cf_bunchers_updated <- iteration_results$cf_bunchers
        c0_updated <- cf_bunchers_updated/bins_bunchers

        # add data to excess_updated_df
        res <- c(as.integer(iteration),bunchers_excess_updated)
        excess_updated_df <- rbind(excess_updated_df,res)

        #  first iteration is 0 (the original count)
        b_diff <- excess_updated_df$bunchers_excess_updated[iteration] -
            excess_updated_df$bunchers_excess_updated[iteration+1]
        iteration <- iteration + 1
    }

    # b_corrected
    b_corrected <- as.numeric(sprintf("%.3f", bunchers_excess_updated/c0_updated))

    # B corrected (excess mass without normalisation)
    B_corrected <- bunchers_excess_updated
    # if we did correction, get last cf (this is the correct one) and assign to cf_graph for graphing
    thedata$cf_density <- iteration_results$cf_density

    # get residuals
    thedata$residuals <- thedata$cf_density - thedata$freq_orig

    # generate output (updated with correction)
    output <- list("data" = thedata, "b_corrected" = b_corrected, "B_corrected" = B_corrected)

    return(output)
}
