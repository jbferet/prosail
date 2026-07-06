#' This function applies noise on a matrix
#'
#' @param refl_lut numeric. reflectance LUT provided as a matrix or dataframe
#' @param noise_level numeric. value of the normal noise proportional to refl_lut
#' to be applied on refl_lut
#' @param noise_type character.
#' - relative: noise proportional to actual value to add noise to
#' - absolute: noise not proportional to actual value to add noise to
#'
#' @return lut_noise numeric. Matrix including data with added noise, same type
#' as refl_lut
#' @export
#' @examples
#' input_prosail <- get_input_prosail(atbd = TRUE, nb_samples = 100)
#' refl_lut <- generate_lut_prosail(input_prosail = input_prosail,
#'                                  spec_prospect = prosail::spec_prospect_full_range,
#'                                  spec_soil = prosail::spec_soil_atbd_v2,
#'                                  spec_atm = prosail::spec_atm)
#' refl_lut_noise <- apply_noise_lut(refl_lut = refl_lut$surf_refl,
#'                                   noise_level = 0.02,
#'                                   noise_type = 'relative')
#'
apply_noise_lut <- function(refl_lut, noise_level = 0, noise_type = 'relative'){
  nb_features <- nrow(refl_lut)
  nb_samples <- ncol(refl_lut)
  if (noise_type == 'relative'){
    lut_noise <- refl_lut + refl_lut*matrix(rnorm(nb_features*nb_samples,0,noise_level),
                                  nrow = nb_features)
  } else if (noise_type == 'absolute'){
    lut_noise <- refl_lut + matrix(rnorm(nb_features*nb_samples,0,noise_level),
                              nrow = nb_features)
  }
  return(lut_noise)
}


#' @rdname prosail-deprecated
#' @export
Apply_Noise_LUT <- function(LUT, NoiseLevel, NoiseType = 'relative'){
  .Deprecated("apply_noise_lut")
  apply_noise_lut(LUT, NoiseLevel, NoiseType)
}


