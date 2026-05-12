#' J4 function for treating (near) conservative scattering
#'
#' @param m numeric. Extinction coefficient for direct (solar or observer) flux
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
#'
jfunc4 <- function(m,t){

  del <- m*t
  out <- 0*del
  out[del>1e-3] <- (1-exp(-del))/(m*(1+exp(-del)))
  out[del<=1e-3] <- 0.5*t*(1.-del*del/12.)
  return(out)
}

#' @rdname prosail-deprecated
#' @export
Jfunc4 <- function(m, t){
  .Deprecated("jfunc4")
  jfunc4(m, t)
}
