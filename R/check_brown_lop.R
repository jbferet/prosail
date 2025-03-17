#' check if brown leaf optical properties are correctly defined
#'
#' @param BrownLOP dataframe. should include wvl, Reflectance & Transmittance
#' @param lambda numeric. spectral bands corresponding to data to simulate
#' @param Input_PROSPECT dataframe. includes all prospect input parameters
#'
#' @return invisible
#' @export

check_brown_lop <- function(BrownLOP, lambda, Input_PROSPECT){
  if  (!'Reflectance' %in%names(BrownLOP) |
       !'Transmittance' %in%names(BrownLOP) |
       !'wvl' %in%names(BrownLOP)){
    message('BrownLOP must include "wvl", "Reflectance" and "Transmittance"')
    stop()
  }
  # spectral domain for brown vegetation matching input optical domain
  if (length(setdiff(lambda, BrownLOP$wvl))>0){
    message('Same spectral domain expected for BrownLOP & SpecPROSPECT')
    stop()
  }
  if (dim(Input_PROSPECT)[1]>1){
    message('BrownLOP defined along with multiple leaf chemical properties')
    message('Only first set of leaf chemical properties will be used to simulate green vegetation')
  }
}
