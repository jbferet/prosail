#' Computes fraction of absorbed photosyntehtically active radiation fAPAR
#'
#' The direct and diffuse light are taken into account as proposed by:
#' Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
#' measurements into total shortwave broadband albedo using a canopy
#' radiative transfer model, Agronomie
#' Es = direct
#' Ed = diffuse
#'
#' @param abs_dir numeric. fraction of direct light absorbed
#' @param abs_hem numeric. fraction of diffuse light absorbed
#' @param tts numeric. Solar zenith angle
#' @param spec_atm_sensor list. direct & diffuse radiation for clear conditions
#' @param par_range numeric. range (in nm) of spectral domain to integrate for
#' computation of fAPAR
#'
#' @return fapar numeric. fAPAR
#' @export
compute_fapar  <- function(abs_dir, abs_hem, tts, spec_atm_sensor,
                           par_range = c(400, 700)){

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
  top <- (abs_dir*par_diro+abs_hem*par_difo)
  par_domain <- which(spec_atm_sensor$lambda >= par_range[1] &
                       spec_atm_sensor$lambda <= par_range[2])
  fapar <- sum(top[par_domain])/sum(par_diro[par_domain]+par_difo[par_domain])
  return(fapar)
}
