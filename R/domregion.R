#' Dominated Region
#'
#' Estimate z (the value of z_vector) that demarcates the upper bound of the dominated region (in notch settings only).

#'
#' @inheritParams bunchit
#' @return  \code{domregion} returns a list with the following objects related to the dominated region (in notch settings only):
#' \item{zD}{The level of z that demarcates the upper bound of the dominated region.}
#' \item{zD_bin}{The value of the bin which zD falls in.}
#' @seealso \code{\link{bunchit}}
#'
#' @examples
#' domregion(zstar = 10000, t0 = 0, t1 = 0.2, binwidth = 50)
#' @export

domregion <- function(zstar,t0,t1, binwidth) {
    zD <- round(zstar*(1-t0)/(1-t1),2)
    zD_bin <- (zD - zstar)/binwidth
    zD <-ceiling(zD)
    zD_bin <- ceiling(zD_bin)

    output <- list("zD" = zD,
                   "zD_bin" = zD_bin)
    return(output)
}


#' @seealso \code{\link{bunchit}}
