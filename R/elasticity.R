#' Elasticity
#'
#' Estimate elasticity from single normalized bunching observation.

#'
#' @param beta normalised excess mass.
#' @inheritParams bunchit
#' @seealso \code{\link{bunchit}},  \code{\link{elasticity_sd}}
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

        if(e_parametric) {
            # define function to suppress iterations output from BBoptim
            hush <- function(code){
                sink("NUL")
                tmp = code
                sink()
                return(tmp)
            }

            # use BB's BBoptim solver to get elasticity
            e <- hush(BBoptim(0.01, notch_equation,
                              t0 = t0, t1 = t1, zstar = zstar, dzstar = Dz))
            e <- e$value
        } else {
            # notch equation (approximation from Kleven's 2018 note, equation 5)
            e <- (1/(2+Dz_over_zstar))*(Dz_over_zstar**2)/(dt/(1-t0))
        }

    }
    return(e)
}
