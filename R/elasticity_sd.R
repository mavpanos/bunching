#' Estimate s.d. of elasticity.

#'
#' @param b_vector vector of normalised excess mass estimates.
#' @inheritParams bunchit
#' @seealso \code{\link{elasticity}} \code{\link{bunchit}}
#' @return  Standard deviation of elasticity estimates.
#' @export

elasticity_sd <- function(b_vector, binwidth, zstar, t0, t1) {
    round(sd(elasticity(b_vector, binwidth, zstar, t0, t1), na.rm = T),3)
}
