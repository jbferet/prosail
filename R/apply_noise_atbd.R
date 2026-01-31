#' This function applies noise defined in S2 ATBD to reflectance
#' look up table
#'
#' @param refl_lut numeric. reflectance look up table
#'
#' @return refl_lut_noise numeric.
#' @export

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
