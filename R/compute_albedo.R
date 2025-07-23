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
#' @param spec_atm_sensor list. direct & diffuse radiation for clear conditions
#' @param par_range numeric. range (in nm) of spectral domain to integrate for
#' computation of albedo
#'
#' @return albedo numeric. albedo
#' @export
compute_albedo  <- function(rsdstar, rddstar, tts, spec_atm_sensor,
                            par_range = c(400, 2400)){

  ############################## #
  ##	direct / diffuse light	##
  ############################## #
  es <- spec_atm_sensor$direct_light
  ed <- spec_atm_sensor$diffuse_light
  rd <- pi/180
  # diffuse radiation (Francois et al., 2002)
  skyl <- 0.847 - 1.61*sin((90-tts)*rd) + 1.04*sin((90-tts)*rd)*sin((90-tts)*rd)
  par_diro <- (1-skyl)*es
  par_difo <- skyl*ed
  top <- (rsdstar*par_diro + rddstar*par_difo)
  albedo_domain <- which(spec_atm_sensor$lambda >= par_range[1] &
                           spec_atm_sensor$lambda <= par_range[2])
  albedo <- sum(top[albedo_domain])/sum(par_diro[albedo_domain] +
                                          par_difo[albedo_domain])
  return(albedo)
}


#' @rdname prosail-deprecated
#' @export
Compute_albedo  <- function(rsdstar, rddstar, tts, SpecATM_Sensor,
                            PAR_range = c(400, 2400)){
  .Deprecated("compute_albedo")
  compute_albedo(rsdstar, rddstar, tts, SpecATM_Sensor, PAR_range)
}
