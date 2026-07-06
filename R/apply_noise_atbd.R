#' This function applies noise defined in S2 ATBD to reflectance
#' look up table
#'
#' @param refl_lut numeric. reflectance look up table provided as a matrix or dataframe
#'
#' @return refl_lut_noise numeric. same type as refl_lut
#' @export
#' @examples
#' input_prosail <- get_input_prosail(atbd = TRUE, nb_samples = 100)
#' refl_lut <- generate_lut_prosail(input_prosail = input_prosail,
#'                                  spec_prospect = prosail::spec_prospect_full_range,
#'                                  spec_soil = prosail::spec_soil_atbd_v2,
#'                                  spec_atm = prosail::spec_atm)
#' refl_lut_noise <- apply_noise_atbd(refl_lut = refl_lut$surf_refl)
#'
apply_noise_atbd <- function(refl_lut){
  ad <- ai <- 0.01
  md <- mi <- 0.02
  wl_b2_b3 <- which(row.names(refl_lut)=='B2' | row.names(refl_lut)=='B3')
  wl_misc <- which(!row.names(refl_lut)=='B2' & !row.names(refl_lut)=='B3')
  refl_lut_noise <- 0*refl_lut

  # add multiplicative noise to B2 and B3
  ad_full <- matrix(rnorm(length(wl_b2_b3)*ncol(refl_lut),0,ad),
                    nrow = length(wl_b2_b3))
  ai_full <- matrix(rnorm(length(wl_b2_b3)*ncol(refl_lut),0,ai),
                    nrow = length(wl_b2_b3))
  md_full <- matrix(rnorm(length(wl_b2_b3)*ncol(refl_lut),0,md),
                    nrow = length(wl_b2_b3))
  MIfull <- matrix(rnorm(length(wl_b2_b3)*ncol(refl_lut),0,mi),
                   nrow = length(wl_b2_b3))
  refl_lut_noise[wl_b2_b3,] <- refl_lut[wl_b2_b3,]*(1+(md_full+MIfull)) +
    ad_full + ai_full

  # add multiplicative noise to other bands
  ad_full <- matrix(rnorm(length(wl_misc)*ncol(refl_lut),0,ad),
                    nrow = length(wl_misc))
  md_full <- matrix(rnorm(length(wl_misc)*ncol(refl_lut),0,md),
                    nrow = length(wl_misc))
  refl_lut_noise[wl_misc,] <- refl_lut[wl_misc,]*(1+md_full) + ad_full
  return(refl_lut_noise)
}
