#' This function converts high resolution spectral info into broader
#' spectral characteristics
#'
#' @param spec_prospect list. optical constants required for PROSPECT
#' @param spec_soil list. dry soil and wet soil, or a unique soil sample if
#'  psoil parameter is not inverted
#' @param spec_atm list. direct and diffuse radiation for clear conditions
#' @param srf list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @return spec_sensor list. list of input specs with sensor resolution
#' @export

prepare_sensor_simulation <- function(spec_prospect, spec_soil, spec_atm, srf){

  # adjust optical constants
  wvl <- spec_prospect$lambda
  # leaf properties
  spec_prospect_sensor <- apply_sensor_characteristics(wvl,spec_prospect,srf)
  # atmospheric properties
  spec_atm_sensor <- apply_sensor_characteristics(wvl,spec_atm,srf)
  # soil properties
  spec_soil_sensor <- apply_sensor_characteristics(wvl,spec_soil,srf)
  spec_sensor <- list('spec_prospect_sensor' = spec_prospect_sensor,
                      'spec_atm_sensor' = spec_atm_sensor,
                      'spec_soil_sensor' = spec_soil_sensor,
                      'band_names' = srf$spectral_bands)
  return(spec_sensor)
}

#' @rdname prosail-deprecated
#' @export
PrepareSensorSimulation <- function(SpecPROSPECT,SpecSOIL,SpecATM,SRF){
  .Deprecated("prepare_sensor_simulation")
  prepare_sensor_simulation(SpecPROSPECT,SpecSOIL,SpecATM,SRF)
}
