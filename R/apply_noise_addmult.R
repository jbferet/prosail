#' This function applies additive and multiplicative noise to reflectance
#' - if additive_noise and multiplicative_noise are defined with a unique value,
#' noise is homogeneous across all spectrum
#' - if additive_noise and multiplicative_noise are the same length as the nb
#' of spectral bands (rows in refl_lut), noise is specific to each spectral band
#'
#' @param refl_lut numeric. reflectance look up table
#' @param additive_noise numeric. additive noise (0 = 0%, 1 = 100%)
#' @param multiplicative_noise numeric. multiplicative noise (0 = 0%, 1 = 100%)
#'
#' @return refl_lut_noise numeric.
#' @export

apply_noise_addmult <- function(refl_lut, additive_noise = 0.01,
                                multiplicative_noise = 0.02){
  nb_wl <- nrow(refl_lut)
  nb_samples <- ncol(refl_lut)
  # add noise to reflectance
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
  refl_lut_noise <- refl_lut*(1+(mult_comp)) + add_comp
  return(refl_lut_noise)
}


#' @rdname prosail-deprecated
#' @export
apply_noise_AddMult <- function(BRF_LUT, AdditiveNoise = 0.01,
                                MultiplicativeNoise = 0.02){
  .Deprecated("apply_noise_addmult")
  apply_noise_addmult(BRF_LUT, AdditiveNoise, MultiplicativeNoise)
}
