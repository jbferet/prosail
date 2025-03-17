#' Computes albedo
#' borrowed from python package PROSAIL developed by Jose Gomez Dans
#'
#' The direct and diffuse light are taken into account as proposed by:
#' Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
#' measurements into total shortwave broadband albedo using a canopy
#' radiative transfer model, Agronomie
#' Es = direct
#' Ed = diffuse
#'
#' @param rsdstar numeric. reflectance from direct light
#' @param rddstar numeric. reflectance from diffuse light
#' @param tts numeric. Solar zenith angle
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param PAR_range numeric. range (in nm) of spectral domain to integrate for
#' computation of albedo
#'
#' @return albedo numeric. albedo
#' @export
compute_albedo  <- function(rsdstar, rddstar, tts, SpecATM_Sensor,
                            PAR_range = c(400, 2400)){

  ############################## #
  ##	direct / diffuse light	##
  ############################## #
  Es <- SpecATM_Sensor$Direct_Light
  Ed <- SpecATM_Sensor$Diffuse_Light
  rd <- pi/180
  # diffuse radiation (Francois et al., 2002)
  skyl <- 0.847 - 1.61*sin((90-tts)*rd) + 1.04*sin((90-tts)*rd)*sin((90-tts)*rd)
  PARdiro <- (1-skyl)*Es
  PARdifo <- skyl*Ed
  top <- (rsdstar*PARdiro + rddstar*PARdifo)
  albedo_domain <- which(SpecATM_Sensor$lambda >= PAR_range[1] &
                           SpecATM_Sensor$lambda <= PAR_range[2])
  albedo <- sum(top[albedo_domain])/sum(PARdiro[albedo_domain] +
                                          PARdifo[albedo_domain])
  return(albedo)
}
