#' Estimate z that demarcates dominated region (in notch settings).

#'
#' @inheritParams bunchit
#' @return  The z level and bin  above zstar where individual is indifferent between there or at zstar. Used for notch
#' @export

domregion <- function(zstar,t0,t1, binwidth) {
    zD <- round(zstar*(1-t0)/(1-t1),0)
    zD_bin <- (zD - zstar)/binwidth
    zD <-ceiling(zD)
    zD_bin <- ceiling(zD_bin)

    output <- list("zD" = zD,
                   "zD_bin" = zD_bin)
    return(output)
}
