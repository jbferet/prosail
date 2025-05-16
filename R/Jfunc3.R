#' J3 function with avoidance of singularity problem
#'
#' @param k numeric. Extinction coefficient for direct (solar or observer) flux
#' @param l numeric.
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
#'
jfunc3 <- function(k,l,t){
  out <- (1.-exp(-(k+l)*t))/(k+l)
  return(out)
}
