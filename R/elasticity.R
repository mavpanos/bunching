#' Estimate elasticity from single normalised bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  Elasticity estimate.
#' @export

elasticity <- function(beta, binwidth, zstar, t0, t1, notch) {
    Dz <- beta*binwidth
    dt <- t1 - t0
    if(notch == F) {
        # kink equation
        e <- -log(1+(Dz/zstar))/log(1-(dt/(1-t0)))
    } else {
        # notch equation
        e <- (Dz/zstar)**2/(dt/(1-t0))
    }
    return(e)
}
