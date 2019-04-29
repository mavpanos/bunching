#' Prepare binned data and model for bunching estimation.
#'
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @export

prepare_data <- function(z_vector, binv, zstar, binwidth, bins_l, bins_r,
                         poly, bins_excl_l,  bins_excl_r, rn, extra_fe) {

    # --------------------------------------------------------------------
    #         generate bin cutoffs based on choice of binv
    # --------------------------------------------------------------------

    # first, get max and min value of running variable, calculate bin number each value belongs to
    zmax <- zstar + (binwidth*bins_r)
    zmin <- zstar - (binwidth*bins_l)
    bins <- seq(zmin,zmax, by = binwidth)

    # generate "thebin" that each observation belongs to, based on binv
    if(binv == "min") {
        thebin <- cut(z_vector, bins, right = F, labels = F)
        thebin <- zmin + binwidth * (thebin - 1)
    } else if(binv == "max") {
        thebin <- cut(z_vector, bins, right = T, labels = F)
        thebin <- zmin + binwidth * thebin
    } else if(binv == "median") {
        thebin <- cut(z_vector, (bins+binwidth/2), right = F, labels = F)
        thebin <- zmin + binwidth * (thebin)
        # in median version, change the maximum bin to NAs since that bin is mechanically only defined over half the binwidth
        thebin[which(thebin == max(thebin, na.rm = T))] <- NA
    }

    # ------------------------------------------------
    #         generate counts per bin
    # ------------------------------------------------
    # turn into dataframe
    thedata <- data.frame("z" = z_vector, "bin" = thebin)
    data_binned <- thedata %>%
        dplyr::group_by(bin) %>%
        dplyr::summarise(freq = n(),
                         z = mean(z, na.rm = T)) %>%
        dplyr::filter(!is.na(bin))

    # add freq_orig column to use after for integration constraint and bootstrap
    data_binned$freq_orig <- data_binned$freq
    data_binned <- as.data.frame(data_binned)


    # checks: if extra kinks do not correspond to a bin, flag it
    if(sum(!is.na(extra_fe) > 0) & (length(extra_fe) != sum(extra_fe %in% data_binned$bin))) {
        print("Warning: extra FE(s) not a valid bin value")
    }


    ##################################################
    #        DATA PREP: bin relative to zstar
    ##################################################

    data_binned$z_rel = (data_binned$bin - zstar)/binwidth

    ##################################################
    #        DATA PREP: dummy for zstar (kink)
    ##################################################

    data_binned$kink <- ifelse(data_binned$bin == zstar,1,0)

    ##################################################
    #  DATA PREP: dummy for extra kinks to control for
    ##################################################

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
    ##################################################
    #   DATA PREP: polynomials
    ##################################################

    polynomial_vector <- c()
    # generate polynomials in z_rel
    for(i in seq(poly)){
        poly_varname <- paste0("poly_",i)
        data_binned[[poly_varname]] <- data_binned$z_rel**i
        polynomial_vector <- c(polynomial_vector, poly_varname)
    }

    ##################################################
    # DATA PREP: dummies for excluded region
    ##################################################

    bins_excluded_all <- c()

    #-----------------
    # below kinkpoint
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
    # above kinkpoint
    #-----------------
    if(bins_excl_r != 0){

        bins_excl_r_vector <- c()

        # otherwise seq(0) gives 1 so it includes a dummy for that one
        for(i in seq(bins_excl_r)) {
            bin_excl_r_varname <- paste0("bin_excl_r_",i)
            data_binned[[bin_excl_r_varname]] <- ifelse(data_binned$z_rel == i,1,0)
            bins_excl_r_vector <- c(bins_excl_r_vector, bin_excl_r_varname)
        }
        # add to bins_excluded_all
        bins_excluded_all <- c(bins_excluded_all,bins_excl_r_vector)
    }

    ##############################################
    # DATA PREP: INDICATOR FOR BUNCHING REGION
    ##############################################

    # make indicator for bunching region
    #  if bins_excluded_all exist (i.e. we didnt set bunching region to just kink)
    # sum bins_excludd_all and kink
    if (length(bins_excluded_all) > 0) {
        data_binned$bunch_region <- rowSums(data_binned[,c("kink",bins_excluded_all)])
    } else {
        data_binned$bunch_region <- data_binned$kink
    }
    ######################################################
    # DATA PREP: INDICATOR FOR BINS ABOVE BUNCHING REGION
    ######################################################
    # create dummy for those right of (above) excluded bunching region
    ul <- zstar + binwidth * bins_excl_r
    data_binned$bin_above_excluded <- ifelse(data_binned$bin > ul,1,0)


    ##############################################
    # DATA PREP: ROUND NUMBER BUNCHING
    ##############################################

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


    ##############################################
    #  EQUATION TO PASS TO LM
    ##############################################

    rhs_vars <- c("kink", extra_fe_vector,
                  polynomial_vector,
                  rn_vector,
                  bins_excluded_all)
    # some may have "" e.g. extra_fe_vector if NA, so drop these
    rhs_vars <- setdiff(rhs_vars, "")
    model_formula <- as.formula(paste0("freq", " ~ ", paste(rhs_vars, collapse = " +")))



    data_forreg <- list("data_binned" = data_binned,
                        "model_formula" = model_formula)
    return(data_forreg)
}
