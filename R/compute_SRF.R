#' Computes spectral response function based on wavelength and fwhm
#' characteristics
#'
#' @param wvl numeric. spectral sampling of the sensor
#' @param fwhm numeric. Full Width Half Maximum for each spectral band
#' @param sensor_name character. sensor name
#'
#' @return SRF list. spectral response, sensor Spectral Bands & Original Bands
#' for which SRF is defined
#' @importFrom stats dnorm
#' @export

compute_SRF <- function(wvl,fwhm, sensor_name = 'Custom'){

  # define full spectral domain in optical domain
  lambda <- prosail::SpecPROSPECT_FullRange$lambda
  Spectral_Response <- matrix(0,ncol = length(wvl),nrow = length(lambda))
  for (i in seq_len(length(wvl))){
    y <- dnorm(lambda,wvl[i],fwhm[i]/2.355)
    Spectral_Response[,i] <- y/max(y)
  }
  Spectral_Response[which(Spectral_Response < 0.001)] <- 0
  SRF <- list('Spectral_Response' = Spectral_Response,
              'Spectral_Bands' = wvl,
              'Original_Bands' = lambda,
              'Central_WL' = wvl,
              'Sensor' = sensor_name)
  return(SRF)
}
