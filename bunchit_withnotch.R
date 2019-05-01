binned_data <- bunching::bin_data(z_vector, binv, zstar, binwidth, bins_l, bins_r)

# -----------------------------
# 2. first pass prep
# -----------------------------

# if it is a notch and we use data-driven method to calculate zu (B = M),
# (re)set bins_excl_r <- 1 for first iteration
if(T > 0 & force_notch == F) {
    bins_excl_r <- 1

    firstpass_prep <- bunching::prep_data_for_reg(binned_data, zstar, binwidth, bins_l, bins_r,
                                                  poly, bins_excl_l, bins_excl_r, rn, extra_fe)
    # start with 1 bin above
    newvar <- paste0("bin_excl_r_", zu_bin)
    firstpass_prep$data_binned[[newvar]]  <- ifelse(firstpass_prep$data_binned$z_rel == zu_bin,1,0)
    firstpass_prep$model_formula <- as.formula(paste(Reduce(paste, deparse(firstpass_prep$model_formula)), newvar, sep = " + "))
    #themodelformula <- as.formula(paste(Reduce(paste, deparse(themodelformula)), newvar, sep = " + "))
    firstpassv2 <- fit_bunching_v2(firstpass_prep$data_binned, firstpass_prep$model_formula)

    # get bunching mass below and missing mass above
    B_below <- firstpassv2$B_zl_zstar
    M_above <- ifelse(!is.na(firstpassv2$M_zstar_zu), firstpassv2$M_zstar_zu, 0)
    # set upper bound on number of iterations: bins above zstar cannot be more than min of
    # 1. available bins
    # 2. remaining DoF
    available_bins <- (max(firstpass_prep$data_binned$bin) - zstar)/binwidth
    DoF_remaining <- dim(firstpass_prep$data_binned)[1] -  dim(firstpassv2$coefficients)[1] - 1
    notch_iterations_bound <- min(available_bins, DoF_remaining)

    zu_bin <- 0
    while ((B_below > M_above) & (zu_bin < notch_iterations_bound)) {
        zu_bin <- zu_bin + 1
        newvar <- paste0("bin_excl_r_", zu_bin)
        firstpass_prep$data_binned[[newvar]]  <- ifelse(firstpass_prep$data_binned$z_rel == zu_bin,1,0)
        firstpass_prep$model_formula <- as.formula(paste(Reduce(paste, deparse(firstpass_prep$model_formula)), newvar, sep = " + "))
        #themodelformula <- as.formula(paste(Reduce(paste, deparse(themodelformula)), newvar, sep = " + "))
        firstpassv2 <- fit_bunching_v2(firstpass_prep$data_binned, firstpass_prep$model_formula)
        #xx2 <- fit_bunching_v2(thedata, themodelformula)
        B_below <- firstpassv2$B_zl_zstar
        M_above <- firstpassv2$M_zstar_zu
    }
}

# AFTER notch calculation, assign results to same object names
# get counterfactual and residuals
counterfactuals_for_graph <- firstpassv2$cf_density
residuals_for_boot <- firstpassv2$residuals
bunchers_initial <- firstpassv2$B_zl_zstar
b_estimate <- firstpassv2$b_estimate_zl_zstar
e_estimate <- NOTCH EQUATION bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1)
model_fit <- firstpassV2$coefficients


} else if (T > 0 & )



if(T == 0) {
firstpass <- bunching::fit_bunching(firstpass_prep$data_binned, firstpass_prep$model_formula)
firstpassv2 <- fit_bunching_v2(firstpass_prep$data_binned, firstpass_prep$model_formula)

} else if (T > 0 & force_notch == F) {
    # Notch analysis




    zu_bin



}

firstpass <- fit_bunching_v2(firstpass_prep$data_binned, firstpass_prep$model_formula)
# get counterfactual and residuals
counterfactuals_for_graph <- firstpass$cf_density
residuals_for_boot <- firstpass$residuals
bunchers_initial <- firstpass$bunchers_excess
b_estimate <- firstpass$b_estimate
e_estimate <- bunching::elasticity(beta = b_estimate, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1)
model_fit <- firstpass$coefficients
# -----------------------------------------
# 3. if no correction needed, do bootstrap
# -----------------------------------------

if(correct == F) {
    boot_results <- bunching::do_bootstrap(firstpass_prep, residuals_for_boot, boot_iterations = n_boot, correction = correct, correction_iterations = iter_max)
    b_sd <- boot_results$b_sd
    b_vector <- boot_results$b_vector
    e_sd <- bunching::elasticity_sd(boot_results$b_vector, binwidth = binwidth, zstar = zstar, t0 = t0, t1 = t1)
    B_for_output <- bunchers_initial # this is bunchers excess. if we dont do integration constraint, this will be output
    B_sd <- boot_results$B_sd
    B_vector <- boot_results$B_vector

}
