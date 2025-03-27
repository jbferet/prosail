#' This function converts high resolution spectral info into broader
#' spectral characteristics
#'
#' @param SpecPROSPECT list. optical constants required for PROSPECT
#' @param SpecSOIL list. dry soil and wet soil, or a unique soil sample if
#'  psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
#' @param SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @return SpecSensor list. list of input specs with sensor resolution
#' @export

prepare_sensor_simulation <- function(SpecPROSPECT,SpecSOIL,SpecATM,SRF){

  # adjust optical constants
  wvl <- SpecPROSPECT$lambda
  # leaf properties
  SpecPROSPECT_Sensor <- apply_sensor_characteristics(wvl,SpecPROSPECT,SRF)
  # atmospheric properties
  SpecATM_Sensor <- apply_sensor_characteristics(wvl,SpecATM,SRF)
  # soil properties
  SpecSOIL_Sensor <- apply_sensor_characteristics(wvl,SpecSOIL,SRF)
  SpecSensor <- list('SpecPROSPECT_Sensor' = SpecPROSPECT_Sensor,
                     'SpecATM_Sensor' = SpecATM_Sensor,
                     'SpecSOIL_Sensor' = SpecSOIL_Sensor,
                     'band_names' = SRF$Spectral_Bands)
  return(SpecSensor)
}
