#' This function adjusts variable values based on co-distribution rules as
#' defined in ATBD co-distributions are all related to lai
#'
#' @param V numeric.
#' @param lai numeric.
#' @param max_lai numeric.
#' @param vmin_0 numeric.
#' @param vmax_0 numeric.
#' @param vmin_lai_max numeric.
#' @param vmax_lai_max numeric.
#'
#' @return Vstar numeric.
#' @export

get_codistributions <- function(V, lai, max_lai, vmin_0, vmax_0,
                                vmin_lai_max, vmax_lai_max){

  vmin_lai <- vmin_0 + (lai*(vmin_lai_max-vmin_0)/max_lai)
  vmax_lai <- vmax_0 + (lai*(vmax_lai_max-vmax_0)/max_lai)
  Vstar <- vmin_lai+((V-vmin_0)*(vmax_lai-vmin_lai)/(vmax_0-vmin_0))
  return(Vstar)
}
