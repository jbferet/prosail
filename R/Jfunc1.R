#' J1 function with avoidance of singularity problem
#'
#' @param k numeric. Extinction coefficient for direct (solar or observer) flux
#' @param l numeric.
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
#'
jfunc1 <- function(k,l,t){
  # J1 function with avoidance of singularity problem
  del <- (k-l)*t
  Jout <- 0*l
  Jout[which(abs(del)>1e-3)] <- (exp(-l[which(abs(del)>1e-3)]*t) -
                                   exp(-k*t))/(k-l[which(abs(del)>1e-3)])
  Jout[which(abs(del)<=1e-3)] <- 0.5*t*(exp(-k*t) +
                                          exp(-l[which(abs(del)<=1e-3)]*t))*
    (1-del[which(abs(del)<=1e-3)]*del[which(abs(del)<=1e-3)]/12)
  return(Jout)
}

#' @rdname prosail-deprecated
#' @export
Jfunc1 <- function(k,l,t){
  .Deprecated("jfunc1")
  jfunc1(k,l,t)
}
