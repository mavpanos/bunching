#' Elasticity s.d.
#' Estimate standard deviation of elasticity.

#'
#' @param b_vector vector of normalised excess mass estimates.
#' @inheritParams bunchit
#' @seealso \code{\link{elasticity}} \code{\link{bunchit}}
#' @return  \code{elasticity_sd} returns the standard deviation of the elasticity estimates.
#' @export

elasticity_sd <- function(b_vector, binwidth, zstar, t0, t1, notch) {
    round(stats::sd(bunching::elasticity(b_vector, binwidth, zstar, t0, t1, notch), na.rm = T),3)
}
