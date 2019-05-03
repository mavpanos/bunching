#' Estimate bootstrapped bunching standard error
#' @param firstpass_prep (binned) data that includes all variables necessary for fitting the model.
#' @param residuals residuals from (first pass) fitted bunchingm model
#' @param boot_iterations number of bootstrapped samples
#' @param correction whether to implement correction for integration constraint.
#' @param correction_iterations maximum iterations for integration constraint correction.
#' @inheritParams bunchit
#'
#' @return b_vector, b_sd, B_vector, B_sd

#' @seealso \code{\link{bunchit}}, \code{\link{prepare_data}}
#' @export


do_bootstrap <- function(firstpass_prep, residuals, boot_iterations, correction,
                         correction_iterations, notch, zD_bin, seed) {
    # set seed if chosen
    if(!is.na(seed)) {
        set.seed(seed)
    }

    # retrieve data and model from firstpass_prep
    data_for_boot <- firstpass_prep$data_binned
    model <- firstpass_prep$model_formula

    # get vector of bootstrapped betas
    boot_results <- sapply(seq(1:boot_iterations), function(i) {
        # adjust frequencies using residual
        data_for_boot$freq_orig <- data_for_boot$freq_orig + sample(residuals, replace = T)
        # make this "freq" so we can pass to fit_bunching which requires "freq ~ ..."
        data_for_boot$freq <- data_for_boot$freq_orig
        # next, re-run first pass on this new series
        booted_firstpass <- bunching::fit_bunching(data_for_boot, model, notch, zD_bin)

        if(correction == F) { # if no need for correction, just take this b
            b_boot <- booted_firstpass$b_estimate
            B_boot <- booted_firstpass$bunchers_excess
            alpha_boot <- booted_firstpass$alpha

        # otherwise, do correction, then take b estimate
        } else if (correction == T) {
            booted_correction <- bunching::do_correction(data_for_boot, booted_firstpass,
                                                         correction_iterations, notch, zD_bin)
            b_boot <- booted_correction$b_corrected
            B_boot <- booted_correction$B_corrected
            alpha_boot <- booted_correction$alpha_corrected
        }
        output <- list("b_boot" = b_boot, "B_boot" = B_boot, "alpha_boot" = alpha_boot)
        return(output)
    })

    b_boot <- unlist(boot_results["b_boot",])
    B_boot <- unlist(boot_results["B_boot",])
    alpha_boot <- unlist(boot_results["alpha_boot",])
    b_sd <- round(sd(b_boot, na.rm = T),3)
    B_sd <- round(sd(B_boot, na.rm = T),3)
    alpha_sd <- round(sd(alpha_boot, na.rm = T),3)

    output <- list("b_vector" = b_boot, "b_sd" = b_sd,
                   "B_vector" = B_boot, "B_sd" = B_sd,
                   "alpha_vector" = alpha_boot, "alpha_sd" = alpha_sd)
    return(output)
}

