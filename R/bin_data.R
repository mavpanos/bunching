#' Bin the raw data
#'
#' Create data frame of binned counts

#'
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  \code{bin_data} returns a data frame with bins and corresponding frequencies (counts).

#' @examples
#' data(bunching_data)
#' binned_data <- bin_data(z_vector = bunching_data$kink, zstar = 10000,
#'                         binwidth = 50, bins_l = 20, bins_r = 20)
#' head(binned_data)


#' @export


bin_data <- function(z_vector, binv = "median", zstar, binwidth, bins_l, bins_r) {

    # --------------------------------------------------------------------
    #         1. generate bin cutoffs based on choice of binv
    # --------------------------------------------------------------------

    # first, get max and min value of z_vector, calculate bin number each value belongs to
    zmax <- zstar + (binwidth*bins_r)
    zmin <- zstar - (binwidth*bins_l)
    bins <- seq(zmin,zmax, by = binwidth)

    # generate "thebin" that each observation belongs to, based on binv
    if(binv == "min") {
        thebin <- cut(z_vector, bins, right = FALSE, labels = FALSE)
        thebin <- zmin + binwidth * (thebin - 1)
    } else if(binv == "max") {
        thebin <- cut(z_vector, bins, right = TRUE, labels = FALSE)
        thebin <- zmin + binwidth * thebin
    } else if(binv == "median") {
        thebin <- cut(z_vector, (bins+binwidth/2), right = FALSE, labels = FALSE)
        thebin <- zmin + binwidth * (thebin)
        # in median version, change the maximum bin to NAs since that bin is
        # mechanically only defined over half the binwidth
        thebin[which(thebin == max(thebin, na.rm = TRUE))] <- NA
    }

    # --------------------------------------------------------------------
    #                   2. generate counts per bin
    # --------------------------------------------------------------------
    # turn into dataframe
    thedata <- data.frame("z" = z_vector, "bin" = thebin)
    thedata <- thedata %>%
        dplyr::group_by(bin) %>%
        dplyr::summarise(freq = n(),
                         z = mean(z, na.rm = TRUE)) %>%
        dplyr::filter(!is.na(bin))

    # add freq_orig column to use after for integration constraint and bootstrap
    thedata$freq_orig <- thedata$freq
    thedata <- as.data.frame(thedata)
    return(thedata)
}
