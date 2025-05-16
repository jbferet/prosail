#' This function applies noise defined in S2 ATBD to BRF LUT
#'
#' @param brf_lut numeric. BRF LUT
#'
#' @return brf_lut_noise numeric.
#' @export

apply_noise_atbd <- function(brf_lut){
  ad <- ai <- 0.01
  md <- mi <- 0.02
  wl_b2_b3 <- which(row.names(brf_lut)=='B2' | row.names(brf_lut)=='B3')
  wl_misc <- which(!row.names(brf_lut)=='B2' & !row.names(brf_lut)=='B3')
  brf_lut_noise <- 0*brf_lut

  # add multiplicative noise to B2 and B3
  ad_full <- matrix(rnorm(length(wl_b2_b3)*ncol(brf_lut),0,ad),
                    nrow = length(wl_b2_b3))
  ai_full <- matrix(rnorm(length(wl_b2_b3)*ncol(brf_lut),0,ai),
                    nrow = length(wl_b2_b3))
  md_full <- matrix(rnorm(length(wl_b2_b3)*ncol(brf_lut),0,md),
                    nrow = length(wl_b2_b3))
  MIfull <- matrix(rnorm(length(wl_b2_b3)*ncol(brf_lut),0,mi),
                   nrow = length(wl_b2_b3))
  brf_lut_noise[wl_b2_b3,] <- brf_lut[wl_b2_b3,]*(1+(md_full+MIfull)) +
    ad_full + ai_full

  # add multiplicative noise to other bands
  ad_full <- matrix(rnorm(length(wl_misc)*ncol(brf_lut),0,ad),
                    nrow = length(wl_misc))
  md_full <- matrix(rnorm(length(wl_misc)*ncol(brf_lut),0,md),
                    nrow = length(wl_misc))
  brf_lut_noise[wl_misc,] <- brf_lut[wl_misc,]*(1+md_full) + ad_full
  return(brf_lut_noise)
}
