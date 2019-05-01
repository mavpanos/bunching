#' Estimate z that demarcates dominated region (in notch settings).

#'
#' @inheritParams bunchit
#' @return  The bin above zstar where individual is indifferent between there or at zstar. Used for notch
#' @export

domregion <- function(zstar,t0,t1, binwidth) {
    zD <- round(zstar*(1-t0)/(1-t1),0)
    ceiling(zD)
    #bin <- (zD - zstar)/binwidth
    #ceiling(bin)
}
