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

compute_srf <- function(wvl,fwhm, sensor_name = 'Custom'){

  # define full spectral domain in optical domain
  lambda <- prosail::spec_prospect_fullrange$lambda
  spectral_response <- matrix(0,ncol = length(wvl),nrow = length(lambda))
  for (i in seq_len(length(wvl))){
    y <- dnorm(lambda,wvl[i],fwhm[i]/2.355)
    spectral_response[,i] <- y/max(y)
  }
  spectral_response[which(spectral_response < 0.001)] <- 0
  srf <- list('Spectral_Response' = spectral_response,
              'Spectral_Bands' = wvl,
              'Original_Bands' = lambda,
              'Central_WL' = wvl,
              'Sensor' = sensor_name)
  return(srf)
}
