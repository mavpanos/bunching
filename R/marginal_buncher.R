#' Marginal Buncher
#'
#' Calculate location (value of z_vector) of marginal buncher.

#'
#' @param beta normalized excess mass.
#' @inheritParams bunchit
#' @param alpha the proportion of individuals in dominated region (in notch setting).
#' @return  \code{marginal_buncher} returns the location of the marginal buncher, i.e. zstar + Dzstar.
#' @seealso \code{\link{bunchit}}
#'
#' @examples
#' marginal_buncher(beta = 2, binwidth = 50, zstar = 10000)
#' @export

marginal_buncher <- function(beta, binwidth, zstar, notch = FALSE, alpha = NULL) {
    if(!notch) {
        # kink specification
        zstar + (beta*binwidth)
    } else {
        # notch specification
        zstar + (beta*binwidth)/(1-alpha)
    }
}

