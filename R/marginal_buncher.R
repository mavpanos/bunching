#' Calculate location (value of running variable) of marginal buncher.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  Returns the location of the marginal buncher, i.e. zstar + Dzstar
#' @export
marginal_buncher <- function(beta, binwidth, zstar) {
    zstar + (beta*binwidth)
}

