#' check if brown leaf optical properties are correctly defined
#'
#' @param brown_lop dataframe. should include wvl, Reflectance & Transmittance
#' @param lambda numeric. spectral bands corresponding to data to simulate
#' @param input_prospect dataframe. includes all prospect input parameters
#'
#' @return invisible
#' @export

check_brown_lop <- function(brown_lop, lambda, input_prospect){
  if (utils::packageVersion("prospect")<'2.0.0'){
    if  (!'Reflectance' %in% names(brown_lop) |
         !'Transmittance' %in% names(brown_lop) |
         !'wvl' %in% names(brown_lop)){
      message('brown_lop must include "wvl", "Reflectance" and "Transmittance"')
      stop()
    }
  } else {
    if  (!'reflectance' %in% names(brown_lop) |
         !'transmittance' %in% names(brown_lop) |
         !'wvl' %in% names(brown_lop)){
      message('brown_lop must include "wvl", "reflectance" and "transmittance"')
      stop()
    }
  }
  # spectral domain for brown vegetation matching input optical domain
  if (length(setdiff(lambda, brown_lop$wvl))>0){
    message('Same spectral domain expected for brown_lop & SpecPROSPECT')
    stop()
  }
  if (dim(input_prospect)[1]>1){
    message('brown_lop defined along with multiple leaf chemical properties')
    message('green vegetation simulated with first set of leaf chemical traits')
  }
}


#' @rdname prosail-deprecated
#' @export
check_BrownLOP <- function(BrownLOP, lambda, Input_PROSPECT){
  .Deprecated("check_brown_lop")
  check_brown_lop(BrownLOP, lambda, Input_PROSPECT)
}

