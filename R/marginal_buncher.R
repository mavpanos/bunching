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
    zstar + (beta*binwidth)
}

