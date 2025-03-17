#' check if spectral sampling is identical between SpecPROSPECT, SpecSOIL,
#' SpecATM
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
#'
#' @return invisible
#' @export
#'
check_spectral_sampling <- function(SpecPROSPECT, SpecSOIL, SpecATM){
  l1 <- SpecPROSPECT$lambda
  l2 <- SpecSOIL$lambda
  l3 <- SpecATM$lambda

  if (!length(l1)==length(l2) | !length(l1)==length(l3) |
      !length(l3)==length(l2)){
    message('Please ensure matching spectral sampling between SpecPROSPECT, SpecSOIL and SpecATM')
    stop()
  } else if (length(unique(l1-l2))>1 | length(unique(l1-l3))>1 |
             length(unique(l3-l2))>1) {
    message('Please ensure matching spectral sampling between SpecPROSPECT, SpecSOIL and SpecATM')
    stop()
  } else if (!unique(l1-l2)==0 | !unique(l1-l3)==0 |
             !unique(l3-l2)==0){
    message('Please ensure matching spectral sampling between SpecPROSPECT, SpecSOIL and SpecATM')
    stop()
  }
  return(invisible())
}
