#' This function applied noise on a matrix
#'
#' @param LUT numeric. Matrix including data to add noise to
#' @param noise_level numeric. value of the normal noise proportional to LUT to
#' apply on LUT
#' @param noise_type character.
#' - relative: noise proportional to actual value to add noise to
#' - absolute: noise not proportional to actual value to add noise to
#'
#' @return LUT_Noise numeric. Matrix including data with added noise
#' @export

apply_noise_lut <- function(LUT, noise_level = 0, noise_type = 'relative'){
  nb_features <- nrow(LUT)
  nb_samples <- ncol(LUT)
  if (noise_type == 'relative'){
    LUT_Noise <- LUT + LUT*matrix(rnorm(nb_features*nb_samples,0,noise_level),
                                  nrow = nb_features)
  } else if (noise_type == 'absolute'){
    LUT_Noise <- LUT + matrix(rnorm(nb_features*nb_samples,0,noise_level),
                              nrow = nb_features)
  }
  return(LUT_Noise)
}
