#' Bootstrap
#'
#' Estimate bunching on bootstrapped samples, using residual-based bootstrapping with replacement.
#' @param firstpass_prep (binned) data that includes all variables necessary for fitting the model.
#' @param residuals residuals from (first pass) fitted bunching model.
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
#'
#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' prepped_data <- prep_data_for_fit(binned_data, zstar = 10000, binwidth = 50,
#'                                   bins_l = 20, bins_r = 20, poly = 4)
#' firstpass <- fit_bunching(prepped_data$data_binned,
#'                           prepped_data$model_formula,
#'                           binwidth = 50)
#' residuals_for_boot <- fit_bunching(prepped_data$data_binned,
#'                                    prepped_data$model_formula,
#'                                    binwidth = 50)$residuals
#' boot_results <- do_bootstrap(zstar = 10000, binwidth = 50,
#'                              firstpass_prep = prepped_data,
#'                              residuals = residuals_for_boot,
#'                              seed = 1)
#' boot_results$b_sd

#' @export


do_bootstrap <- function(zstar, binwidth, firstpass_prep, residuals, n_boot = 100,
                         correct = TRUE, correct_iter_max = 200, notch = FALSE, zD_bin = NA, seed = NA) {
    # set seed if chosen
    if(!is.na(seed)) {
        set.seed(seed)
    }

    # retrieve data and model from firstpass_prep
    data_for_boot <- firstpass_prep$data_binned
    model <- firstpass_prep$model_formula

    # get vector of bootstrapped betas
    boot_results <- sapply(seq(1:n_boot), function(i) {
        # adjust frequencies using residual
        data_for_boot$freq_orig <- data_for_boot$freq_orig + sample(residuals, replace = TRUE)
        # make this "freq" so we can pass to fit_bunching which requires "freq ~ ..."
        data_for_boot$freq <- data_for_boot$freq_orig
        # next, re-run first pass on this new series
        booted_firstpass <- bunching::fit_bunching(data_for_boot, model, binwidth, notch, zD_bin)

        # if no need for integration correction, just take this b
        if(correct == FALSE) {
            b_boot <- booted_firstpass$b_estimate
            B_boot <- booted_firstpass$bunchers_excess
            alpha_boot <- booted_firstpass$alpha
            mbuncher_boot <- bunching::marginal_buncher(beta = b_boot, binwidth = binwidth, zstar = zstar)

        # otherwise, do correction, then take b estimate
        } else if (correct == TRUE) {
            booted_correction <- bunching::do_correction(zstar, binwidth, data_for_boot, booted_firstpass,
                                                         correct_iter_max, notch, zD_bin)
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
    b_sd <- stats::sd(b_boot, na.rm = TRUE)
    B_sd <- stats::sd(B_boot, na.rm = TRUE)
    alpha_sd <- stats::sd(alpha_boot, na.rm = TRUE)
    mbuncher_sd <- stats::sd(mbuncher_boot, na.rm = TRUE)

    output <- list("b_vector" = b_boot, "b_sd" = b_sd,
                   "B_vector" = B_boot, "B_sd" = B_sd,
                   "marginal_buncher_vector" = mbuncher_boot, "marginal_buncher_sd" = mbuncher_sd,
                   "alpha_vector" = alpha_boot, "alpha_sd" = alpha_sd)
    return(output)
}

