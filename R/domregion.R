#' Estimate z (the value of the running variable) that demarcates the dominated region (in notch settings only).

#'
#' @inheritParams bunchit
#' @return  The level of z, and the bin it is in (above zstar), that demarcates the dominated region (in notch settings only).
#' @export

domregion <- function(zstar,t0,t1, binwidth) {
    zD <- round(zstar*(1-t0)/(1-t1),2)
    zD_bin <- (zD - zstar)/binwidth
    zD <-ceiling(zD)
    zD_bin <- ceiling(zD_bin)

    output <- list("zD" = zD,
                   "zD_bin" = zD_bin)
    return(output)
}
