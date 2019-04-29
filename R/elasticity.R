#' Estimate elasticity from single bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  Elasticity estimate.
#' @export

elasticity <- function(beta, binwidth, zstar, t0, t1) {
    Dz <- beta*binwidth
    dt <- t1-t0
    e <-  (Dz/zstar)/((t1-t0)/(1-t0))
    return(e)
}
