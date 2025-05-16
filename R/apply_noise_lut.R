#' This function applied noise on a matrix
#'
#' @param lut numeric. Matrix including data to add noise to
#' @param noise_level numeric. value of the normal noise proportional to lut to
#' apply on lut
#' @param noise_type character.
#' - relative: noise proportional to actual value to add noise to
#' - absolute: noise not proportional to actual value to add noise to
#'
#' @return lut_noise numeric. Matrix including data with added noise
#' @export

apply_noise_lut <- function(lut, noise_level = 0, noise_type = 'relative'){
  nb_features <- nrow(lut)
  nb_samples <- ncol(lut)
  if (noise_type == 'relative'){
    lut_noise <- lut + lut*matrix(rnorm(nb_features*nb_samples,0,noise_level),
                                  nrow = nb_features)
  } else if (noise_type == 'absolute'){
    lut_noise <- lut + matrix(rnorm(nb_features*nb_samples,0,noise_level),
                              nrow = nb_features)
  }
  return(lut_noise)
}
