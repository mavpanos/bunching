#' Indifference Level in Notch Setting
#'
#' Calculate location (value of running variable) where marginal buncher is indifferent between staying there and moving to notch.

#'
#' @param notch_marginal_buncher location of marginal buncher in notch, i.e. zI for the one initially at zstar + Dzstar.
#' @inheritParams bunchit
#' @return  \code{notch_indifference} returns the indifference location of the marginal buncher, i.e. zI.
#' @seealso \code{\link{bunchit}}
#' @export


notch_indifference <- function(notch_marginal_buncher, t0, t1, elasticity) {
    notch_marginal_buncher * ((1-t1)/(1-t0))^elasticity
}
