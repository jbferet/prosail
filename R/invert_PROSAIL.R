#' Performs PROSAIL inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param brf_mes numeric. measured BRF to be used during inversion
#' @param initialization list. Initial guess of Parms2Estimate for the inversion
#' @param lower_bound list. Lower bound
#' @param upper_bound list. Upper bound
#' @param SpecPROSPECT_Sensor  list. Includes optical constants
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param SpecSOIL_Sensor list. includes either dry soil and wet soil,
#' or a unique soil sample if the psoil parameter is not inverted
#' @param TypeLidf numeric. Type of leaf inclination distribution function
#' @param parm_set  list. Parameters set to a fixed value by user
#' @param merit_function  character. name of the merit function
#' with given criterion to minimize (default = RMSE)
#' @param prior_info list. prior mean, sd and weight of parameters defined as xprior
#' modulate its importance
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom pracma fmincon
#' @export

invert_PROSAIL  <- function(brf_mes, initialization = NULL, lower_bound,
                            upper_bound, SpecPROSPECT_Sensor, SpecATM_Sensor,
                            SpecSOIL_Sensor, TypeLidf, parm_set,
                            merit_function = 'merit_RMSE_PROSAIL',
                            prior_info = NULL){

  # define parameters included in inversion or set
  ParmInv <- which_parms_to_invert(initialization = initialization,
                                    lower_bound = lower_bound,
                                    upper_bound = upper_bound,
                                    parm_set = parm_set)
  Parms2Prior <- NULL
  if (!is.null(prior_info)){
    prior_info_sort <- which_parm_prior(prior_mean = prior_info$mean,
                                prior_sd = prior_info$sd)
    prior_info$mean <- prior_info_sort$mean
    prior_info$sd <- prior_info_sort$sd
    Parms2Prior <- prior_info_sort$Parms2Prior
  }
  ParmInv$InVar[ParmInv$Parms2Set] <- ParmInv$parm_set

  xinit <- as.numeric(ParmInv$initialization)
  lb <- as.numeric(ParmInv$lower_bound)
  names(lb) <- names(ParmInv$lower_bound)
  ub <- as.numeric(ParmInv$upper_bound)
  names(ub) <- names(ParmInv$upper_bound)
  resInv <- fmincon(x0 = xinit, fn = merit_function, gr = NULL,
                    parms_xinit = names(ParmInv$initialization),
                    brf_mes = brf_mes,
                    SpecPROSPECT_Sensor = SpecPROSPECT_Sensor,
                    SpecSOIL_Sensor = SpecSOIL_Sensor,
                    SpecATM_Sensor = SpecATM_Sensor,
                    Parms2Estimate = ParmInv$Parms2Estimate,
                    input_prosail = ParmInv$InVar, TypeLidf = TypeLidf,
                    prior_info = prior_info, Parms2Prior = Parms2Prior,
                    method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                    lb = lb, ub = ub, hin = NULL, heq = NULL, tol = 1e-14,
                    maxfeval = 2000, maxiter = 2000)
  ParmInv$InVar[ParmInv$Parms2Estimate] <- resInv$par
  return(ParmInv$InVar)
}
