#' This function applied noise on a matrix
#'
#' @param LUT numeric. Matrix including data to add noise to
#' @param NoiseLevel numeric. value of the normal noise proportional to LUT to
#' apply on LUT
#' @param NoiseType character.
#' - relative: noise proportional to actual value to add noise to
#' - absolute: noise not proportional to actual value to add noise to
#'
#' @return LUT_Noise numeric. Matrix including data with added noise
#' @export

apply_noise_lut <- function(LUT, NoiseLevel, NoiseType = 'relative'){
  nbFeatures <- nrow(LUT)
  nbSamples <- ncol(LUT)
  if (NoiseType == 'relative'){
    LUT_Noise <- LUT + LUT*matrix(rnorm(nbFeatures*nbSamples,0,NoiseLevel),
                                  nrow = nbFeatures)
  } else if (NoiseType == 'absolute'){
    LUT_Noise <- LUT + matrix(rnorm(nbFeatures*nbSamples,0,NoiseLevel),
                              nrow = nbFeatures)
  }
  return(LUT_Noise)
}
