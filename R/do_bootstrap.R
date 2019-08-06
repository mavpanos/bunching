#' Bootstrap
#'
#' Estimate bunching on bootstrapped samples, using residual-based bootstrapping with replacement.
#' @param firstpass_prep (binned) data that includes all variables necessary for fitting the model.
#' @param residuals residuals from (first pass) fitted bunching model.
#' @param boot_iterations number of bootstrapped samples.
#' @param correction whether to implement correction for integration constraint.
#' @param correction_iterations maximum iterations for integration constraint correction.
#' @inheritParams bunchit
#' @inheritParams fit_bunching
#'
#' @return \code{do_bootstrap} returns a list with the following bootstrapped estimates:
#' \item{b_vector}{A vector with the bootstrapped normalized excess mass estimates.}
#' \item{b_sd}{The standard deviation of the bootstrapped b_vector.}
#' \item{B_vector}{A vector with the bootstrapped excess mass estimates (not normalized).}
#' \item{B_sd}{The standard deviation of the bootstrapped B_vector.}
#' \item{marginal_buncher_vector}{A vector with the bootstrapped estimates of the location (z value) of the marginal buncher.}
#' \item{marginal_buncher_sd}{The standard deviation of the bootstrapped marginal_buncher_vector.}
#' \item{alpha_vector}{A vector with the bootstrapped estimates of the fraction of bunchers in the dominated region (only in notch case).}
#' \item{alpha_vector_sd}{The standard deviation of the bootstrapped alpha_vector.}


#' @seealso \code{\link{bunchit}}, \code{\link{prep_data_for_fit}}
#' @export


do_bootstrap <- function(zstar, binwidth, firstpass_prep, residuals, boot_iterations,
                         correction, correction_iterations, notch, zD_bin, seed) {
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

        # if no need for integration correction, just take this b
        if(correction == F) {
            b_boot <- booted_firstpass$b_estimate
            B_boot <- booted_firstpass$bunchers_excess
            alpha_boot <- booted_firstpass$alpha
            mbuncher_boot <- bunching::marginal_buncher(beta = b_boot, binwidth = binwidth, zstar = zstar)

        # otherwise, do correction, then take b estimate
        } else if (correction == T) {
            booted_correction <- bunching::do_correction(zstar, binwidth, data_for_boot, booted_firstpass,
                                                         correction_iterations, notch, zD_bin)
            b_boot <- booted_correction$b_corrected
            B_boot <- booted_correction$B_corrected
            alpha_boot <- booted_correction$alpha_corrected
            mbuncher_boot <- bunching::marginal_buncher(beta = b_boot, binwidth = binwidth, zstar = zstar)
        }
        output <- list("b_boot" = b_boot, "B_boot" = B_boot,
                       "alpha_boot" = alpha_boot, "mbuncher_boot" = mbuncher_boot)
        return(output)
    })

    b_boot <- unlist(boot_results["b_boot",])
    B_boot <- unlist(boot_results["B_boot",])
    alpha_boot <- unlist(boot_results["alpha_boot",])
    mbuncher_boot <- unlist(boot_results["mbuncher_boot",])
    b_sd <- stats::sd(b_boot, na.rm = T)
    B_sd <- stats::sd(B_boot, na.rm = T)
    alpha_sd <- stats::sd(alpha_boot, na.rm = T)
    mbuncher_sd <- stats::sd(mbuncher_boot, na.rm = T)

    output <- list("b_vector" = b_boot, "b_sd" = b_sd,
                   "B_vector" = B_boot, "B_sd" = B_sd,
                   "marginal_buncher_vector" = mbuncher_boot, "marginal_buncher_sd" = mbuncher_sd,
                   "alpha_vector" = alpha_boot, "alpha_sd" = alpha_sd)
    return(output)
}

