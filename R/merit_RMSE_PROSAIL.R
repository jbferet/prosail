#' Merit function for PROSAIL inversion
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param parms_xinit character. name of parameters corresponding to xinit
#' @param brfMES  numeric. measured BRF
#' @param SpecPROSPECT_Sensor list. Includes optical constants for PROSPECT
#' @param SpecSOIL_Sensor list. Includes soil reflectance (either 2 references
#' for optimization of psoil, or a unique spectrun)
#' @param SpecATM_Sensor list. Includes diffuse and direct light
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param Parms2Estimate  numeric. rank of variables to be inverted
#' to be estimated through inversion
#' @param Parm2Set  numeric. rank of variables to be set out of inversion
#' @param ParmSet  list. value of variables to be set out of inversion
#' @param InVar dataframe. full set of PROSAIL input variables
#' @param TypeLidf  numeric. type of leaf inclination distribution function
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior SD of parameters defined as xprior
#' @param Parms2Prior numeric. rank of parameters to be used with prior info
#' @param WeightPrior numeric. Weight to be applied on prior information
#' to modulate its importance
#'
#' @return fc estimates of the parameters
#' @export
merit_RMSE_PROSAIL <- function(xinit, parms_xinit, brfMES, SpecPROSPECT_Sensor,
                               SpecSOIL_Sensor, SpecATM_Sensor, Parms2Estimate,
                               Parm2Set, ParmSet, InVar, TypeLidf,
                               PriorInfoMean = NULL, PriorInfoSD = NULL,
                               Parms2Prior = NULL, WeightPrior = 0.01){

  xinit[xinit<0] <- 0
  InVar[Parms2Estimate] <- xinit[match(parms_xinit, Parms2Estimate)]
  xprior <- InVar[Parms2Prior]
  rsoil <- InVar$psoil*SpecSOIL_Sensor$Dry_Soil +
    (1-InVar$psoil)*SpecSOIL_Sensor$Wet_Soil
  # call PROSAIL to get reflectance from 4 fluxes
  Ref <- PRO4SAIL(Spec_Sensor = SpecPROSPECT_Sensor, Input_PROSPECT = InVar,
                  TypeLidf = TypeLidf, LIDFa = InVar$LIDFa, LIDFb = InVar$LIDFb,
                  lai = InVar$lai, q = InVar$q, tts = InVar$tts,
                  tto = InVar$tto, psi = InVar$psi, rsoil = rsoil)
  # Computes BRF based on outputs from PROSAIL and sun position
  brfMOD <- compute_BRF(rdot = Ref$rdot, rsot = Ref$rsot,
                        tts = InVar$tts, SpecATM_Sensor = SpecATM_Sensor)
  # compute cost
  fc <- cost_function_RMSE_PROSAIL(brfMES = brfMES, brfMOD$BRF, xprior,
                             PriorInfoMean = PriorInfoMean,
                             PriorInfoSD = PriorInfoSD,
                             WeightPrior = WeightPrior)
  return(fc)
}
