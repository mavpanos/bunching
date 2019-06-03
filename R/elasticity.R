#' Elasticity
#'
#' Estimate elasticity from single normalized bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}},  \code{\link{elasticity_sd}}
#' @return  \code{elasticity} returns the estimated elasticity. In the kink case it uses the standard log-form. In the notch case, it is the reduced-form approximation in equation (5) of Kleven's 2018 note "Calculating Reduced-Form Elasticities Using Notches".
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
