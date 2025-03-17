#' Performs PROSAIL inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param brfMES numeric. measured BRF to be used during inversion
#' @param InitialGuess list. Initial guess of Parms2Estimate for the inversion
#' @param LowerBound list. Lower bound
#' @param UpperBound list. Upper bound
#' @param SpecPROSPECT_Sensor  list. Includes optical constants
#' refractive index, specific absorption coefs and corresponding spectral bands
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param SpecSOIL_Sensor list. includes either dry soil and wet soil,
#' or a unique soil sample if the psoil parameter is not inverted
#' @param TypeLidf numeric. Type of leaf inclination distribution function
#' @param ParmSet  list. Parameters set to a fixed value by user
#' @param MeritFunction  character. name of the merit function
#' with given criterion to minimize (default = RMSE)
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior SD of parameters defined as xprior
#' @param WeightPrior numeric. Weight to be applied on prior information to
#' modulate its importance
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom pracma fmincon
#' @export

invert_PROSAIL  <- function(brfMES, InitialGuess = NULL, LowerBound, UpperBound,
                            SpecPROSPECT_Sensor, SpecATM_Sensor,
                            SpecSOIL_Sensor, TypeLidf, ParmSet,
                            MeritFunction = 'merit_RMSE_PROSAIL',
                            PriorInfoMean = NULL, PriorInfoSD = NULL,
                            WeightPrior = 0.01){

  # define parameters included in inversion or set
  ParmInv <- which_parms_to_invert(InitialGuess = InitialGuess,
                                    LowerBound = LowerBound,
                                    UpperBound = UpperBound,
                                    ParmSet = ParmSet)
  Parms2Prior <- NULL
  if (!is.null(PriorInfoMean) & !is.null(PriorInfoSD)){
    PriorInfo <- which_parm_prior(PriorInfoMean = PriorInfoMean,
                                PriorInfoSD = PriorInfoSD)
    PriorInfoMean <- PriorInfo$PriorInfoMean
    PriorInfoSD <- PriorInfo$PriorInfoSD
    Parms2Prior <- PriorInfo$Parms2Prior
  }
  ParmInv$InVar[ParmInv$Parms2Set] <- ParmInv$ParmSet

  xinit <- as.numeric(ParmInv$InitialGuess)
  lb <- as.numeric(ParmInv$LowerBound)
  names(lb) <- names(ParmInv$LowerBound)
  ub <- as.numeric(ParmInv$UpperBound)
  names(ub) <- names(ParmInv$UpperBound)
  resInv <- fmincon(x0 = xinit, fn = MeritFunction, gr = NULL,
                    parms_xinit = names(ParmInv$InitialGuess), brfMES = brfMES,
                    SpecPROSPECT_Sensor = SpecPROSPECT_Sensor,
                    SpecSOIL_Sensor = SpecSOIL_Sensor,
                    SpecATM_Sensor = SpecATM_Sensor,
                    Parms2Estimate = ParmInv$Parms2Estimate,
                    Parm2Set = ParmInv$Parms2Set, ParmSet = ParmInv$ParmSet,
                    InVar = ParmInv$InVar, TypeLidf = TypeLidf,
                    PriorInfoMean = PriorInfoMean, PriorInfoSD = PriorInfoSD,
                    Parms2Prior = Parms2Prior, WeightPrior = WeightPrior,
                    method = "SQP",A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                    lb = lb, ub = ub, hin = NULL, heq = NULL, tol = 1e-14,
                    maxfeval = 2000, maxiter = 2000)
  ParmInv$InVar[ParmInv$Parms2Estimate] <- resInv$par
  return(ParmInv$InVar)
}
