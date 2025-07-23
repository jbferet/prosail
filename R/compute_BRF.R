#' Computes BRF based on outputs from PROSAIL & SZA
#'
#' The direct and diffuse light are taken into account as proposed by:
#' Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
#' measurements into total shortwave broadband albedo using a canopy
#' radiative transfer model, Agronomie
#' Es = direct
#' Ed = diffuse
#'
#' @param rdot numeric. Hemispherical-directional refl factor viewing direction
#' @param rsot numeric. Bi-directional reflectance factor
#' @param tts numeric. Solar zenith angle
#' @param spec_atm_sensor list. direct & diffuse radiation for clear conditions
#' @param skyl numeric. values for skyl
#' @return BRF numeric. Bidirectional reflectance factor
#' @export

compute_brf  <- function(rdot, rsot, tts, spec_atm_sensor, skyl = NULL){

  ############################## #
  ##	direct / diffuse light	##
  ############################## #
  es <- spec_atm_sensor$direct_light
  ed <- spec_atm_sensor$diffuse_light
  rd <- pi/180
  # diffuse radiation (Francois et al., 2002)
  if (is.null(skyl))
    skyl <- 0.847-1.61*sin((90-tts)*rd) + 1.04*sin((90-tts)*rd)*sin((90-tts)*rd)
  par_diro <- (1-skyl)*es
  par_difo <- skyl*ed
  brf <- (rdot*par_difo+rsot*par_diro)/(par_diro+par_difo)
  return(data.frame('BRF' = brf))
}


#' @rdname prosail-deprecated
#' @export
Compute_BRF  <- function(rdot, rsot, tts, SpecATM_Sensor, skyl = NULL){
  .Deprecated("compute_brf")
  compute_brf(rdot, rsot, tts, SpecATM_Sensor, skyl)
}
