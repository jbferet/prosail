#' This function applies additive and multiplicative noise to BRF data
#' - if additive_noise and multiplicative_noise are defined with a unique value,
#' noise is homogeneous across all spectrum
#' - if additive_noise and multiplicative_noise are the same length as the nb
#' of spectral bands (rows in brf_lut), noise is specific to each spectral band
#'
#' @param brf_lut numeric. BRF LUT
#' @param additive_noise numeric. additive noise (0 = 0%, 1 = 100%)
#' @param multiplicative_noise numeric. multiplicative noise (0 = 0%, 1 = 100%)
#'
#' @return brf_lut_noise numeric.
#' @export

apply_noise_addmult <- function(brf_lut, additive_noise = 0.01,
                                multiplicative_noise = 0.02){
  nb_wl <- nrow(brf_lut)
  nb_samples <- ncol(brf_lut)
  # add noise to BRF
  if (length(additive_noise)==1){
    add_comp <- matrix(rnorm(nb_wl*nb_samples,0,additive_noise),
                       nrow = nb_wl)
  } else if ((length(additive_noise)==nb_wl)){
    add_comp <- matrix(data = 0, nrow = nb_wl, ncol = nb_samples)
    for (i in seq_len(nb_wl))
      add_comp[i,] <- matrix(data = rnorm(nb_samples, mean = 0,
                                          sd = additive_noise[i]),
                             ncol = nb_samples)
  }
  if (length(multiplicative_noise)==1){
    mult_comp <- matrix(rnorm(nb_wl*nb_samples,0,multiplicative_noise),
                        nrow = nb_wl)
  } else if ((length(multiplicative_noise)==nb_wl)){
    mult_comp <- matrix(data = 0, nrow = nb_wl, ncol = nb_samples)
    for (i in seq_len(nb_wl))
      mult_comp[i,] <- matrix(data = rnorm(nb_samples, mean = 0,
                                           sd = multiplicative_noise[i]),
                              ncol = nb_samples)
  }
  brf_lut_noise <- brf_lut*(1+(mult_comp)) + add_comp
  return(brf_lut_noise)
}
