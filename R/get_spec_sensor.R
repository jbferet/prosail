#' Computes spectral sensor characteristics for optical constants required for
#' PROSPECT soil parameters as well as direct and diffuse radiation for clear
#' conditions
#' @param SensorName character. name of the sensor. Either already defined, or
#' to be defined based on other inputs
#' @param SpectralProps list. Sensor spectral characteristics including
#' wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character. Path where to get or save SRF
#'
#' @return SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#' @export

get_spec_sensor <- function(SensorName = 'MyCustomSensor',
                            SpectralProps = NULL,
                            Path_SensorResponse = NULL){
  SRF <- get_radiometry(SensorName,SpectralProps = SpectralProps,
                       Path_SensorResponse = Path_SensorResponse)
  if (is.null(SRF)) stop()
  # adjust optical constants from 1 nm sampling to S2 spectral sampling
  SpecSensor <- prepare_sensor_simulation(prosail::SpecPROSPECT_FullRange,
                                        prosail::SpecSOIL,
                                        prosail::SpecATM,SRF)
  return(SpecSensor)
}
