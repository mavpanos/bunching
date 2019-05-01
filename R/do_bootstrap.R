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


do_bootstrap <- function(firstpass_prep, residuals, boot_iterations, correction, correction_iterations, notch) {
    data_for_boot <- firstpass_prep$data_binned
    model <- firstpass_prep$model_formula

    b_vector <- sapply(seq(1:boot_iterations), function(i) {
        # adjust frequencies using residual
        data_for_boot$freq_orig <- data_for_boot$freq_orig + sample(residuals, replace = T)
        # make this "freq" so we can pass to run_reg which requires "freq ~ ..."
        data_for_boot$freq <- data_for_boot$freq_orig
        # next, re-run first pass on this new series, get residuals out of this
        booted_firstpass <- bunching::fit_bunching(data_for_boot, model, notch)

        if(correction == F) { # if no need for correction, just take this b
            b_boot <- booted_firstpass$b_estimate
            B_boot <- booted_firstpass$bunchers_excess
        } else if (correction == T) {
            booted_correction <- bunching::do_correction(thedata = data_for_boot, firstpass_results = booted_firstpass,
                                               max_iterations = correction_iterations, notch = notch)
            b_boot <- booted_correction$b_corrected
            B_boot <- booted_correction$B_corrected
        }
        b_results <- list("b_boot" = b_boot, "B_boot" = B_boot)
        return(b_results)
    })

    b_boot <- unlist(b_vector["b_boot",])
    B_boot <- unlist(b_vector["B_boot",])
    b_sd <- round(sd(b_boot, na.rm = T),3)
    B_sd <- round(sd(B_boot, na.rm = T),3)

    output <- list("b_vector" = b_boot, "b_sd" = b_sd,
                   "B_vector" = B_boot, "B_sd" = B_sd)
    return(output)
}

