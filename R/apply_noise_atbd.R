#' This function applies noise defined in S2 ATBD to BRF LUT
#'
#' @param BRF_LUT numeric. BRF LUT
#'
#' @return BRF_LUT_Noise numeric.
#' @export

apply_noise_atbd <- function(BRF_LUT){
  AD <- AI <- 0.01
  MD <- MI <- 0.02
  WL_B2_B3 <- which(row.names(BRF_LUT)=='B2' | row.names(BRF_LUT)=='B3')
  WL_misc <- which(!row.names(BRF_LUT)=='B2' & !row.names(BRF_LUT)=='B3')
  BRF_LUT_Noise <- 0*BRF_LUT

  # add multiplicative noise to B2 and B3
  ADfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,AD),
                   nrow = length(WL_B2_B3))
  AIfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,AI),
                   nrow = length(WL_B2_B3))
  MDfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,MD),
                   nrow = length(WL_B2_B3))
  MIfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,MI),
                   nrow = length(WL_B2_B3))
  BRF_LUT_Noise[WL_B2_B3,] <- BRF_LUT[WL_B2_B3,]*(1+(MDfull+MIfull)) +
    ADfull + AIfull

  # add multiplicative noise to other bands
  ADfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,AD),
                   nrow = length(WL_misc))
  AIfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,AI),
                   nrow = length(WL_misc))
  MDfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,MD),
                   nrow = length(WL_misc))
  MIfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,MI),
                   nrow = length(WL_misc))
  BRF_LUT_Noise[WL_misc,] <- BRF_LUT[WL_misc,]*(1+MDfull) + ADfull
  return(BRF_LUT_Noise)
}
