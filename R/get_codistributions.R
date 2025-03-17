#' This function adjusts variable values based on co-distribution rules as
#' defined in ATBD co-distributions are all related to LAI
#'
#' @param V numeric.
#' @param LAI numeric.
#' @param MaxLAI numeric.
#' @param Vmin0 numeric.
#' @param Vmax0 numeric.
#' @param VminLAImax numeric.
#' @param VmaxLAImax numeric.
#'
#' @return Vstar numeric.
#' @export

get_codistributions <- function(V, LAI, MaxLAI, Vmin0, Vmax0,
                                VminLAImax, VmaxLAImax){

  VminLAI <- Vmin0 + (LAI*(VminLAImax-Vmin0)/MaxLAI)
  VmaxLAI <- Vmax0 + (LAI*(VmaxLAImax-Vmax0)/MaxLAI)
  Vstar <- VminLAI+((V-Vmin0)*(VmaxLAI-VminLAI)/(Vmax0-Vmin0))
  return(Vstar)
}
