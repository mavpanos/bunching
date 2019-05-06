#' Estimate elasticity from single normalised bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  Elasticity estimate.
#' @export

elasticity <- function(beta, binwidth, zstar, t0, t1, notch) {
    Dz <- beta*binwidth
    Dz_over_zstar <- Dz/zstar
    dt <- t1 - t0

    if(notch == F) {
        # kink equation
        e <- -log(1+Dz_over_zstar)/log(1-(dt/(1-t0)))
    } else {
        # notch equation (approximation from Kleven's 2018 note, equation 5)
        e <- (1/(2+Dz_over_zstar))*(Dz_over_zstar**2)/(dt/(1-t0))
    }
    return(e)
}
