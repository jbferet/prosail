#' Merit function for PROSAIL inversion
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param parms_xinit character. name of parameters corresponding to xinit
#' @param brf_mes  numeric. measured BRF
#' @param SpecPROSPECT_Sensor list. Includes optical constants for PROSPECT
#' @param SpecSOIL_Sensor list. Includes soil reflectance (either 2 references
#' for optimization of psoil, or a unique spectrun)
#' @param SpecATM_Sensor list. Includes diffuse and direct light
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param Parms2Estimate  numeric. rank of variables to be inverted
#' to be estimated through inversion
#' @param input_prosail dataframe. full set of PROSAIL input variables
#' @param TypeLidf  numeric. type of leaf inclination distribution function
#' @param prior_info list. prior mean, sd & weight of parms defined as xprior
#' @param Parms2Prior numeric. rank of parameters to be used with prior info
#' to modulate its importance
#'
#' @return fc estimates of the parameters
#' @export
merit_RMSE_PROSAIL <- function(xinit, parms_xinit, brf_mes, SpecPROSPECT_Sensor,
                               SpecSOIL_Sensor, SpecATM_Sensor, Parms2Estimate,
                               input_prosail, TypeLidf,
                               prior_info = NULL, Parms2Prior = NULL){

  xinit[xinit<0] <- 0
  input_prosail[Parms2Estimate] <- xinit[match(parms_xinit, Parms2Estimate)]
  xprior <- input_prosail[Parms2Prior]
  rsoil <- input_prosail$psoil*SpecSOIL_Sensor$Dry_Soil +
    (1-input_prosail$psoil)*SpecSOIL_Sensor$Wet_Soil
  # call PROSAIL to get reflectance from 4 fluxes
  Ref <- PRO4SAIL(Spec_Sensor = SpecPROSPECT_Sensor,
                  input_prospect = input_prosail,
                  TypeLidf = TypeLidf, LIDFa = input_prosail$LIDFa,
                  LIDFb = input_prosail$LIDFb,
                  lai = input_prosail$lai, q = input_prosail$q,
                  tts = input_prosail$tts, tto = input_prosail$tto,
                  psi = input_prosail$psi, rsoil = rsoil)
  # Computes BRF based on outputs from PROSAIL and sun position
  brf_mod <- compute_BRF(rdot = Ref$rdot, rsot = Ref$rsot,
                         tts = input_prosail$tts,
                         SpecATM_Sensor = SpecATM_Sensor)
  # compute cost
  fc <- cost_function_RMSE_PROSAIL(brf_mes = brf_mes, brf_mod$BRF, xprior,
                                   prior_info = prior_info)
  return(fc)
}
