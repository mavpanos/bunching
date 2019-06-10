#' Marginal Buncher
#'
#' Calculate location (value of running variable) of marginal buncher.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @return  \code{marginal_buncher} returns the location of the marginal buncher, i.e. zstar + Dzstar.
#' @seealso \code{\link{bunchit}}
#' @export

marginal_buncher <- function(beta, binwidth, zstar) {
    if(beta < 0) {
        stop("Bunching estimate is negative. Are you sure this is a kink/notch? \n
             Check your choices of zstar, bins_excl_l and bins_excl_r")
    }
    zstar + (beta*binwidth)
}

