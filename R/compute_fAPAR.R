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
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param PAR_range numeric. range (in nm) of spectral domain to integrate for
#' computation of fAPAR
#'
#' @return fAPAR numeric. fAPAR
#' @export
compute_fAPAR  <- function(abs_dir, abs_hem, tts, SpecATM_Sensor,
                           PAR_range = c(400, 700)){

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
  top <- (abs_dir*PARdiro+abs_hem*PARdifo)
  PARdomain <- which(SpecATM_Sensor$lambda >= PAR_range[1] &
                       SpecATM_Sensor$lambda <= PAR_range[2])
  fAPAR <- sum(top[PARdomain])/sum(PARdiro[PARdomain]+PARdifo[PARdomain])
  return(fAPAR)
}
