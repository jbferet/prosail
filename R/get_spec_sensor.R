#' Computes spectral sensor characteristics for optical constants required for
#' PROSPECT soil parameters as well as direct and diffuse radiation for clear
#' conditions
#' @param sensor_name character. name of the sensor. Either already defined, or
#' to be defined based on other inputs
#' @param spectral_properties list. Sensor spectral characteristics including
#' wl = central wavelength and fwhm = corresponding fwhm
#' @param srf_path character. Path where to get or save SRF
#'
#' @return srf list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#' @export

get_spec_sensor <- function(sensor_name = 'user_defined',
                            spectral_properties = NULL,
                            srf_path = NULL){
  srf <- get_spectral_response_function(sensor_name = sensor_name,
                        spectral_properties = spectral_properties,
                        srf_path = srf_path)
  if (is.null(srf)) stop()
  # adjust optical constants from 1 nm sampling to S2 spectral sampling
  spec_sensor <- prepare_sensor_simulation(spec_prospect = prosail::spec_prospect_fullrange,
                                           spec_soil = prosail::spec_soil_ossl,
                                           spec_atm = prosail::spec_atm,
                                           srf = srf)
  return(spec_sensor)
}
