#' Bin the raw data
#'
#' Create dataframe of binned counts

#'
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  \code{bin_data} returns a data frame with two columns:
#' \item{bin}{The bin value.}
#' \item{freq}{The count per bin.}
#' @export


bin_data <- function(z_vector, binv, zstar, binwidth, bins_l, bins_r) {

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
    thedata <- thedata %>%
        dplyr::group_by(bin) %>%
        dplyr::summarise(freq = n(),
                         z = mean(z, na.rm = T)) %>%
        dplyr::filter(!is.na(bin))

    # add freq_orig column to use after for integration constraint and bootstrap
    thedata$freq_orig <- thedata$freq
    thedata <- as.data.frame(thedata)
    return(thedata)
}
