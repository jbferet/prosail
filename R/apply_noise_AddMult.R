#' This function applies additive and multiplicative noise to BRF data
#' - if AdditiveNoise and MultiplicativeNoise are defined with a unique value,
#' noise is homogeneous across all spectrum
#' - if AdditiveNoise and MultiplicativeNoise are the same length as the number
#' of spectral bands (rows in BRF_LUT), noise is specific to each spectral band
#'
#' @param BRF_LUT numeric. BRF LUT
#' @param AdditiveNoise numeric. additive noise (0 = 0%, 1 = 100%)
#' @param MultiplicativeNoise numeric. multiplicative noise (0 = 0%, 1 = 100%)
#'
#' @return BRF_LUT_Noise numeric.
#' @export

apply_noise_AddMult <- function(BRF_LUT, AdditiveNoise = 0.01,
                                MultiplicativeNoise = 0.02){
  nbWL <- nrow(BRF_LUT)
  nbSamples <- ncol(BRF_LUT)
  # add noise to BRF
  if (length(AdditiveNoise)==1){
    AddComp <- matrix(rnorm(nbWL*nbSamples,0,AdditiveNoise),
                      nrow = nbWL)
  } else if ((length(AdditiveNoise)==nbWL)){
    AddComp <- matrix(data = 0, nrow = nbWL, ncol = nbSamples)
    for (i in seq_len(nbWL))
      AddComp[i,] <- matrix(data = rnorm(nbSamples,
                                         mean = 0,
                                         sd = AdditiveNoise[i]),
                            ncol = nbSamples)
  }
  if (length(MultiplicativeNoise)==1){
    MultComp <- matrix(rnorm(nbWL*nbSamples,0,MultiplicativeNoise), nrow = nbWL)
  } else if ((length(MultiplicativeNoise)==nbWL)){
    MultComp <- matrix(data = 0, nrow = nbWL, ncol = nbSamples)
    for (i in seq_len(nbWL))
      MultComp[i,] <- matrix(data = rnorm(nbSamples,
                                          mean = 0,
                                          sd = MultiplicativeNoise[i]),
                             ncol = nbSamples)
  }
  BRF_LUT_Noise <- BRF_LUT*(1+(MultComp)) + AddComp
  return(BRF_LUT_Noise)
}
