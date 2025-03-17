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
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param skyl numeric. values for skyl
#' @return BRF numeric. Bidirectional reflectance factor
#' @export
compute_BRF  <- function(rdot, rsot, tts, SpecATM_Sensor, skyl = NULL){

  ############################## #
  ##	direct / diffuse light	##
  ############################## #
  Es <- SpecATM_Sensor$Direct_Light
  Ed <- SpecATM_Sensor$Diffuse_Light
  rd <- pi/180
  # diffuse radiation (Francois et al., 2002)
  if (is.null(skyl))
    skyl <- 0.847 - 1.61*sin((90-tts)*rd) +
    1.04*sin((90-tts)*rd)*sin((90-tts)*rd)
  PARdiro <- (1-skyl)*Es
  PARdifo <- skyl*Ed
  BRF <- (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo)
  return(data.frame('BRF' = BRF))
}
