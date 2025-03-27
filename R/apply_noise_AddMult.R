#' This function applies additive and multiplicative noise to BRF data
#' - if additive_noise and multiplicative_noise are defined with a unique value,
#' noise is homogeneous across all spectrum
#' - if additive_noise and multiplicative_noise are the same length as the number
#' of spectral bands (rows in BRF_LUT), noise is specific to each spectral band
#'
#' @param BRF_LUT numeric. BRF LUT
#' @param additive_noise numeric. additive noise (0 = 0%, 1 = 100%)
#' @param multiplicative_noise numeric. multiplicative noise (0 = 0%, 1 = 100%)
#'
#' @return BRF_LUT_Noise numeric.
#' @export

apply_noise_AddMult <- function(BRF_LUT, additive_noise = 0.01,
                                multiplicative_noise = 0.02){
  nb_wl <- nrow(BRF_LUT)
  nb_samples <- ncol(BRF_LUT)
  # add noise to BRF
  if (length(additive_noise)==1){
    AddComp <- matrix(rnorm(nb_wl*nb_samples,0,additive_noise),
                      nrow = nb_wl)
  } else if ((length(additive_noise)==nb_wl)){
    AddComp <- matrix(data = 0, nrow = nb_wl, ncol = nb_samples)
    for (i in seq_len(nb_wl))
      AddComp[i,] <- matrix(data = rnorm(nb_samples, mean = 0,
                                         sd = additive_noise[i]),
                            ncol = nb_samples)
  }
  if (length(multiplicative_noise)==1){
    MultComp <- matrix(rnorm(nb_wl*nb_samples,0,multiplicative_noise),
                       nrow = nb_wl)
  } else if ((length(multiplicative_noise)==nb_wl)){
    MultComp <- matrix(data = 0, nrow = nb_wl, ncol = nb_samples)
    for (i in seq_len(nb_wl))
      MultComp[i,] <- matrix(data = rnorm(nb_samples, mean = 0,
                                          sd = multiplicative_noise[i]),
                             ncol = nb_samples)
  }
  BRF_LUT_Noise <- BRF_LUT*(1+(MultComp)) + AddComp
  return(BRF_LUT_Noise)
}
