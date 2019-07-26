#' Elasticity
#'
#' Estimate elasticity from single normalized bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}}
#' @return  \code{elasticity} returns the estimated elasticity. By default, this is based on the reduced-form approximation. To use the parametric equivalent, set e_parametric to TRUE.
#' @export

elasticity <- function(beta, binwidth, zstar, t0, t1, notch, e_parametric) {

    # define some quantities to simplify equations
    Dz <- beta*binwidth
    Dz_over_zstar <- Dz/zstar
    dt <- t1 - t0

    # kinks elasticity, parametric and reduced form
    if(notch == F) {
        if(e_parametric) {
            e <- -log(1+Dz_over_zstar)/log(1-(dt/(1-t0)))
        } else {
            e <- Dz_over_zstar/(dt/(1-t0))
        }



    } else {
        # notch elasticity, parametric and reduced form

        #  first, calculate reduced-form for both cases
        # we use this if parametri does not converge
        e <- (1/(2+Dz_over_zstar))*(Dz_over_zstar**2)/(dt/(1-t0))


        if(e_parametric) {
            # silence iterations' output of BBoptim by sinking messages to tmp file f
            f = file()
            sink(file = f)

            # use BB's BBoptim solver to get elasticity
            # catch any errors, and suppress warnings
            lower_bound <- .00001
            upper_bound <- 10
            estimate <- tryCatch({
                suppressWarnings(BB::BBoptim(0.01, notch_equation,
                                             t0 = t0, t1 = t1, zstar = zstar, dzstar = Dz,
                                             lower = lower_bound, upper = upper_bound))
            },
            error=function(error_message) {
                return(error_message)
            })
            # switch off silencing, close file
            sink()
            close(f)

            # if convergence field exists in output, algorithm ran (maybe converged or not)
            warning_message <- "The elasticity estimate based on the parametric version for notches has no solution, reduced-form estimate is returned instead."

            if(!"convergence" %in% names(estimate)) {
                # if it doesn't exist, could not run estimation.
                # return warning and keep reduced-form
                warning(warning_message)
            } else {
                if(estimate$convergence != 0) {
                    # if estimated but did not converge (non-0 flag) same
                    warning(warning_message)
                } else {
                    e <- estimate$par
                    # if e hit the upper bound allowed, flag it
                    if(e == upper_bound) {
                        warning("The elasticity estimate based on the parametric version for notches hit the upper bound of possible solution values. Interpet with caution!")
                    }
                    if(e == lower_bound) {
                        warning("The elasticity estimate based on the parametric version for notches hit the lower bound of possible solution values. Interpet with caution!")
                    }
                }
            }
        }




    }


    return(e)
}
