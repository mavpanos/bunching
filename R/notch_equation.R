notch_equation <- function(e, t0, t1, zstar, dzstar) {
    # define some intermediate variable to simplify equation
    one_over_one_plus_dz_over_z <- 1/(1+ (dzstar/zstar))
    delta_t_over_t <- (t1-t0)/(1-t0)

    y <- one_over_one_plus_dz_over_z -
        (1/(1+(1/e)) * (one_over_one_plus_dz_over_z^(1+(1/e)))) -
        ((1/(1+e)) * ((1-delta_t_over_t)^(1+e)))
    y
}

