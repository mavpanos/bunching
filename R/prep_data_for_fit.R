#' Data Preparation
#'
#' Prepare binned data and model for bunching estimation.
#'
#'@param data_binned dataframe of counts per bin
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return \code{data_binned} returns a list with the following:
#' \item{data_binned}{The binned data with the extra columns necessary for model fitting, such as indicators for bunching region, fixed effects, etc.}
#' \item{model_formula}{The formula used for model fitting.}
#'
#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' prepped_data <- prep_data_for_fit(binned_data, zstar = 10000, binwidth = 50,
#'                                   bins_l = 20, bins_r = 20, poly = 4,
#'                                   bins_excl_l = 2, bins_excl_r = 3,
#'                                   rn = c(250,500), extra_fe = 10200)
#' head(prepped_data$data_binned)
#' prepped_data$model_formula

#' @export

prep_data_for_fit <- function(data_binned, zstar, binwidth, bins_l, bins_r, poly = 9,
                              bins_excl_l = 0,  bins_excl_r = 0, rn = NA, extra_fe = NA,
                              correct_above_zu = FALSE) {

    # --------------------------------------------
    # bin relative to zstar
    # --------------------------------------------
    data_binned$z_rel = (data_binned$bin - zstar)/binwidth

    # --------------------------------------------
    # dummy for zstar
    # --------------------------------------------
    data_binned$zstar <- ifelse(data_binned$bin == zstar,1,0)

    # --------------------------------------------
    # dummy for extra kinks to control for
    # --------------------------------------------
    if(sum(!is.na(extra_fe)) > 0) { # default is extra_fe = NA
        # store names of extra_fe_vectors
        extra_fe_vector <- c()
        for(i in extra_fe){
            extra_fe_varname <- paste0("extra_fe_",i)
            data_binned[[extra_fe_varname]] <- ifelse(data_binned$bin == i,1,0)
            extra_fe_vector <- c(extra_fe_vector, extra_fe_varname)
        }
    } else {
        extra_fe_vector <- ""
    }

    # --------------------------------------------
    # polynomials
    # --------------------------------------------
    polynomial_vector <- c()
    # generate polynomials in z_rel
    for(i in seq(poly)){
        poly_varname <- paste0("poly_",i)
        data_binned[[poly_varname]] <- data_binned$z_rel**i
        polynomial_vector <- c(polynomial_vector, poly_varname)
    }

    # --------------------------------------------
    # dummies for excluded region
    # --------------------------------------------
    bins_excluded_all <- c()

    #-----------------
    # below zstar
    #-----------------
    if(bins_excl_l != 0){
        bins_excl_l_vector <- c()
        for(i in seq(bins_excl_l)) {
            bins_excl_l_varname <- paste0("bin_excl_l_",i)
            data_binned[[bins_excl_l_varname]] <- ifelse(data_binned$z_rel == -i,1,0)
            bins_excl_l_vector <- c(bins_excl_l_vector, bins_excl_l_varname)
        }
        # add all bins_excl_l_vector to bins_excluded_all
        bins_excluded_all <- c(bins_excluded_all,bins_excl_l_vector)
    }

    #-----------------
    # above zstar
    #-----------------
    if(bins_excl_r != 0){

        bins_excl_r_vector <- c()

        for(i in seq(bins_excl_r)) {
            bin_excl_r_varname <- paste0("bin_excl_r_",i)
            data_binned[[bin_excl_r_varname]] <- ifelse(data_binned$z_rel == i,1,0)
            bins_excl_r_vector <- c(bins_excl_r_vector, bin_excl_r_varname)
        }
        # add to bins_excluded_all
        bins_excluded_all <- c(bins_excluded_all,bins_excl_r_vector)
    }

    # --------------------------------------------
    # indicator for bunching region
    # --------------------------------------------
    # if bins_excluded_all exist (i.e. we didnt set bunching region to just zstar)
    # sum bins_excludd_all and zstar
    if (length(bins_excluded_all) > 0) {
        data_binned$bunch_region <- rowSums(data_binned[,c("zstar",bins_excluded_all)])
    } else {
        data_binned$bunch_region <- data_binned$zstar
    }

    # --------------------------------------------
    # indicator for bins above bunching region
    # --------------------------------------------
    if(correct_above_zu) {
        # create dummy for those right of (above) excluded bunching region/zstar
        ul <- zstar + binwidth * bins_excl_r
        data_binned$bin_above_excluded <- ifelse(data_binned$bin > ul,1,0)
    } else {
        # else make this all bins above zstar
        data_binned$bin_above_excluded <- ifelse(data_binned$bin > zstar,1,0)
    }

    # --------------------------------------------
    # round number bunching indicators
    # --------------------------------------------
    if(sum(!is.na(rn)) > 0) {
        # make dummies for round numbers and store variable names
        # first, sort them (useful for fixing colinearity below)
        rn <- rn[order(rn)]

        rn_vector <- c()
        for(i in rn){
            rn_dummy_name <- paste0("rn_",i)
            data_binned[[rn_dummy_name]] <- ifelse(data_binned$bin %% i == 0,1,0)
            rn_vector <- c(rn_vector, rn_dummy_name)
        }
        # if two round numbers selected, and one is multiple of the other, replace smaller with 0
        if(length(rn) == 2 & (rn[2] %% rn[1] == 0)) {
            data_binned[[rn_vector[1]]] <- ifelse(data_binned[[rn_vector[1]]] == 1 & data_binned[[rn_vector[2]]] == 1,
                                                  0, data_binned[[rn_vector[1]]])
        }
    } else {
        rn_vector <- ""
    }

    # --------------------------------------------
    #  equation to pass to lm
    # --------------------------------------------

    rhs_vars <- c("zstar", extra_fe_vector,
                  polynomial_vector,
                  rn_vector,
                  bins_excluded_all)
    # some may have "" e.g. extra_fe_vector if NA, so drop these
    rhs_vars <- setdiff(rhs_vars, "")
    model_formula <- stats::as.formula(paste0("freq", " ~ ", paste(rhs_vars, collapse = " +")))



    data_forreg <- list("data_binned" = data_binned,
                        "model_formula" = model_formula)
    return(data_forreg)
}
