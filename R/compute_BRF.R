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
#' @param spec_atm_sensor list. direct and diffuse radiation for clear conditions
#' @param skyl numeric. values for skyl
#' @return BRF numeric. Bidirectional reflectance factor
#' @export
compute_brf <- function(rdot, rsot, tts, spec_atm_sensor, skyl = NULL) {
  ##############################
  ## 	direct / diffuse light	##
  ##############################
  e_s <- spec_atm_sensor$Direct_Light
  e_d <- spec_atm_sensor$Diffuse_Light
  rd <- pi / 180
  # diffuse radiation (Francois et al., 2002)
  if (is.null(skyl)) {
    skyl <- 0.847 - 1.61 * sin((90 - tts) * rd) + 1.04 * sin((90 - tts) * rd) * sin((90 - tts) * rd)
  }
  par_dir_o <- (1 - skyl) * e_s
  par_dif_o <- skyl * e_d
  brf <- (rdot * par_dif_o + rsot * par_dir_o) / (par_dir_o + par_dif_o)
  return(data.frame("BRF" = brf))
}

#' @rdname prosail-deprecated
#' @export
Compute_BRF <- function(rdot, rsot, tts, SpecATM_Sensor, skyl = NULL) {
  .Deprecated("compute_brf")
  compute_brf(rdot, rsot, tts, SpecATM_Sensor, skyl)
}
