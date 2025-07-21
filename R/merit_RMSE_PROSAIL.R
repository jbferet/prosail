#' Merit function for PROSAIL inversion
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param parms_xinit character. name of parameters corresponding to xinit
#' @param brf_mes  numeric. measured BRF
#' @param spec_prospect_sensor list. Includes optical constants for PROSPECT
#' @param spec_soil_sensor list. Includes soil reflectance (either 2 references
#' for optimization of psoil, or a unique spectrun)
#' @param spec_atm_sensor list. Includes diffuse and direct light
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param parms_to_estimate  numeric. rank of variables to be inverted
#' to be estimated through inversion
#' @param input_prosail dataframe. full set of PROSAIL input variables
#' @param type_lidf  numeric. type of leaf inclination distribution function
#' @param prior_info list. prior mean, sd & weight of parms defined as xprior
#' @param parms_to_prior numeric. rank of parameters to be used with prior info
#' to modulate its importance
#'
#' @return fc estimates of the parameters
#' @export
merit_rmse_prosail <- function(xinit, parms_xinit, brf_mes,
                               spec_prospect_sensor, spec_soil_sensor,
                               spec_atm_sensor, parms_to_estimate,
                               input_prosail, type_lidf, prior_info = NULL,
                               parms_to_prior = NULL){

  xinit[xinit<0] <- 0
  input_prosail[parms_to_estimate] <- xinit[match(parms_xinit,
                                                  parms_to_estimate)]
  xprior <- input_prosail[parms_to_prior]
  rsoil <- input_prosail$psoil*spec_soil_sensor$max_refl +
    (1-input_prosail$psoil)*spec_soil_sensor$min_refl
  # call PROSAIL to get reflectance from 4 fluxes
  refl <- prosail(spec_sensor = spec_prospect_sensor,
                  input_prospect = input_prosail,
                  type_lidf = type_lidf, lidf_a = input_prosail$lidf_a,
                  lidf_b = input_prosail$lidf_b,
                  lai = input_prosail$lai, hotspot = input_prosail$hotspot,
                  tts = input_prosail$tts, tto = input_prosail$tto,
                  psi = input_prosail$psi, rsoil = rsoil)
  # Computes BRF based on outputs from PROSAIL and sun position
  brf_mod <- compute_brf(rdot = refl$rdot, rsot = refl$rsot,
                         tts = input_prosail$tts,
                         spec_atm_sensor = spec_atm_sensor)
  # compute cost
  fc <- cost_function_rmse_prosail(brf_mes = brf_mes, brf_mod$BRF, xprior,
                                   prior_info = prior_info)
  return(fc)
}
