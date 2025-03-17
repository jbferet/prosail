#' dcum function
#' @param a numeric. controls the average leaf slope
#' @param b numeric. controls the distribution's bimodality
#' @param t numeric. angle
#' @return f
#' @export
dcum <- function(a,b,t){
  rd <- pi/180
  if (a>=1){
    f <- 1-cos(rd*t)
  } else {
    eps <- 1e-8
    delx <- 1
    x <- 2*rd*t
    p <- x
    while (delx >= eps){
      y <- a*sin(x)+.5*b*sin(2.*x)
      dx <- .5*(y-x+p)
      x <- x+dx
      delx <- abs(dx)
    }
    f <- (2.*y+p)/pi
  }
  return(f)
}
