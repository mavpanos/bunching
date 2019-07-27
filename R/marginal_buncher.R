#' Marginal Buncher
#'
#' Calculate location (value of running variable) of marginal buncher.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @param alpha the proportion of individuals in dominated region (in notch setting).
#' @return  \code{marginal_buncher} returns the location of the marginal buncher, i.e. zstar + Dzstar.
#' @seealso \code{\link{bunchit}}
#' @export

marginal_buncher <- function(beta, binwidth, zstar, notch = F, alpha = NULL) {
    if(!notch) {
        # kink specification
        zstar + (beta*binwidth)
    } else {
        # notch specification
    zstar + (beta*binwidth)/(1-alpha)
    }
}

