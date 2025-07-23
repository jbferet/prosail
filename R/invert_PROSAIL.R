#' Performs PROSAIL inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param brf_mes numeric. measured BRF to be used during inversion
#' @param initialization list. Initial guess of parms_to_estimate
#' @param lower_bound list. Lower bound
#' @param upper_bound list. Upper bound
#' @param spec_prospect_sensor  list. Includes optical constants
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param spec_atm_sensor list. direct & diffuse radiation for clear conditions
#' @param spec_soil_sensor list. includes either dry soil and wet soil,
#' or a unique soil sample if the psoil parameter is not inverted
#' @param type_lidf numeric. Type of leaf inclination distribution function
#' @param parm_set  list. Parameters set to a fixed value by user
#' @param merit_function  character. name of the merit function
#' with given criterion to minimize (default = RMSE)
#' @param prior_info list. prior mean, sd and weight of parameters defined as
#' xprior modulate its importance
#'
#' @return OutPROSPECT estimated values corresponding to parms_to_estimate
#' @importFrom pracma fmincon
#' @export

invert_prosail  <- function(brf_mes, initialization = NULL, lower_bound,
                            upper_bound, spec_prospect_sensor, spec_atm_sensor,
                            spec_soil_sensor, type_lidf, parm_set,
                            merit_function = 'merit_rmse_prosail',
                            prior_info = NULL){

  # define parameters included in inversion or set
  parms_to_invert <- which_parms_to_invert(initialization = initialization,
                                           lower_bound = lower_bound,
                                           upper_bound = upper_bound,
                                           parm_set = parm_set)
  parms_to_prior <- NULL
  if (!is.null(prior_info)){
    prior_info_sort <- which_parm_prior(prior_mean = prior_info$mean,
                                        prior_sd = prior_info$sd)
    prior_info$mean <- prior_info_sort$mean
    prior_info$sd <- prior_info_sort$sd
    parms_to_prior <- prior_info_sort$parms_to_prior
  }
  parms_to_invert$input_prosail[parms_to_invert$parms_to_set] <- parms_to_invert$parm_set

  xinit <- as.numeric(parms_to_invert$initialization)
  lb <- as.numeric(parms_to_invert$lower_bound)
  names(lb) <- names(parms_to_invert$lower_bound)
  ub <- as.numeric(parms_to_invert$upper_bound)
  names(ub) <- names(parms_to_invert$upper_bound)
  res_inv <- fmincon(x0 = xinit, fn = merit_function, gr = NULL,
                     parms_xinit = names(parms_to_invert$initialization),
                     brf_mes = brf_mes,
                     spec_prospect_sensor = spec_prospect_sensor,
                     spec_soil_sensor = spec_soil_sensor,
                     spec_atm_sensor = spec_atm_sensor,
                     parms_to_estimate = parms_to_invert$parms_to_estimate,
                     input_prosail = parms_to_invert$input_prosail,
                     type_lidf = type_lidf,
                     prior_info = prior_info, parms_to_prior = parms_to_prior,
                     method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                     lb = lb, ub = ub, hin = NULL, heq = NULL, tol = 1e-14,
                     maxfeval = 2000, maxiter = 2000)
  parms_to_invert$input_prosail[parms_to_invert$parms_to_estimate] <- res_inv$par
  return(parms_to_invert$input_prosail)
}


#' @rdname prosail-deprecated
#' @export
Invert_PROSAIL  <- function(brfMES, InitialGuess = NULL, LowerBound, UpperBound,
                            SpecPROSPECT_Sensor, SpecATM_Sensor, SpecSOIL_Sensor,
                            TypeLidf, ParmSet, MeritFunction = 'Merit_RMSE_PROSAIL',
                            PriorInfoMean = NULL, PriorInfoSD = NULL,
                            WeightPrior = 0.01){
  .Deprecated("invert_prosail")
  prior_info <- list('mean' = PriorInfoMean,
                     'SD' = PriorInfoSD,
                     'weight_prior' = WeightPrior)
  invert_prosail(brfMES, InitialGuess, LowerBound, UpperBound,
                 SpecPROSPECT_Sensor, SpecATM_Sensor, SpecSOIL_Sensor,
                 TypeLidf, ParmSet, MeritFunction, prior_info)
}
