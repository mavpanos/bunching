#' Notch Equation
#'
#' Defines indifference condition based on parametric utility function in notch setting.

#'
#' @inheritParams bunchit
#' @param e elasticity.
#' @param dzstar The distance of the marginal buncher above zstar.
#' @return  \code{util_diff} returns the difference in utility between zstar and z_I in notch setting.
#' @seealso \code{\link{bunchit}}
#' @seealso \code{\link{elasticity}}
#' @export

notch_equation <- function(e, t0, t1, zstar, dzstar) {
    # define some intermediate variable to simplify equation
    one_over_one_plus_dz_over_z <- 1/(1+ (dzstar/zstar))
    delta_t_over_t <- (t1-t0)/(1-t0)

    util_diff <- one_over_one_plus_dz_over_z -
        (1/(1+(1/e)) * (one_over_one_plus_dz_over_z^(1+(1/e)))) -
        ((1/(1+e)) * ((1-delta_t_over_t)^(1+e)))
    util_diff
}

