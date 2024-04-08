# ============================================================================== =
# prosail
# Lib_PROSAIL_Inversion.R
# ============================================================================== =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================== =
# This Library includes functions dedicated to PROSAIL inversion using iterative
# optimization
# ============================================================================== =

#' Performs PROSAIL inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param brfMES numeric. measured BRF to be used during inversion
#' @param InitialGuess list. Initial guess of Parms2Estimate for the inversion
#' @param LowerBound list. Lower bound
#' @param UpperBound list. Upper bound
#' @param SpecPROSPECT_Sensor  list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param SpecSOIL_Sensor list. includes either dry soil and wet soil,
#' or a unique soil sample if the psoil parameter is not inverted
#' @param TypeLidf numeric. Type of leaf inclination distribution function
#' @param ParmSet  list. Parameters set to a fixed value by user
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard deviation of parameters defined as xprior
#' @param WeightPrior numeric. Weight to be applied on prior information to
#' modulate its importance
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom pracma fmincon
#' @export

Invert_PROSAIL  <- function(brfMES, InitialGuess = NULL, LowerBound, UpperBound,
                            SpecPROSPECT_Sensor, SpecATM_Sensor, SpecSOIL_Sensor,
                            TypeLidf, ParmSet, MeritFunction = 'Merit_RMSE_PROSAIL',
                            PriorInfoMean = NULL, PriorInfoSD = NULL,
                            WeightPrior = 0.01){

  # define parameters included in inversion or set
  ParmInv <- WhichParameters2Invert(InitialGuess = InitialGuess,
                                    LowerBound = LowerBound,
                                    UpperBound = UpperBound,
                                    ParmSet = ParmSet)
  Parms2Prior <- NULL
  if (!is.null(PriorInfoMean) & !is.null(PriorInfoSD)){
    PriorInfo <- WhichParmPrior(PriorInfoMean = PriorInfoMean,
                                PriorInfoSD = PriorInfoSD)
    PriorInfoMean <- PriorInfo$PriorInfoMean
    PriorInfoSD <- PriorInfo$PriorInfoSD
    Parms2Prior <- PriorInfo$Parms2Prior
  }
  ParmInv$InVar[ParmInv$Parms2Set] <- ParmInv$ParmSet

  # update init value and lower/upper boundaries for inversion based on Vars2Estimate
  # xinit <- as.vector(ParmInv$InitialGuess,mode='numeric')
  # lb <- as.vector(ParmInv$LowerBound,mode='numeric')
  # ub <- as.vector(ParmInv$UpperBound,mode='numeric')
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

#' Merit function for PROSAIL inversion
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param parms_xinit character. name of parameters corresponding to xinit
#' @param brfMES  numeric. measured BRF
#' @param SpecPROSPECT_Sensor list. Includes optical constants for PROSPECT
#' @param SpecSOIL_Sensor list. Includes soil reflectance (either 2 references
#' for optimization of psoil, or a unique spectrun)
#' @param SpecATM_Sensor list. Includes diffuse and direct light
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Parms2Estimate  numeric. rank of variables to be inverted
#' to be estimated through inversion
#' @param Parm2Set  numeric. rank of variables to be set out of inversion
#' @param ParmSet  list. value of variables to be set out of inversion
#' @param InVar dataframe. full set of PROSAIL input variables
#' @param TypeLidf  numeric. type of leaf inclination distribution function
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard deviation of parameters defined as xprior
#' @param Parms2Prior numeric. rank of the parameters which should be used with prior information
#' @param WeightPrior numeric. Weight to be applied on prior information to modulate its importance
#'
#' @return fc estimates of the parameters
#' @export
Merit_RMSE_PROSAIL <- function(xinit, parms_xinit, brfMES, SpecPROSPECT_Sensor,
                               SpecSOIL_Sensor, SpecATM_Sensor, Parms2Estimate,
                               Parm2Set, ParmSet, InVar, TypeLidf,
                               PriorInfoMean = NULL, PriorInfoSD = NULL,
                               Parms2Prior = NULL, WeightPrior = 0.01){

  xinit[xinit<0] <- 0
  InVar[Parms2Estimate] <- xinit[match(parms_xinit, Parms2Estimate)]
  xprior <- InVar[Parms2Prior]
  rsoil <- InVar$psoil*SpecSOIL_Sensor$Dry_Soil+(1-InVar$psoil)*SpecSOIL_Sensor$Wet_Soil
  # call PROSAIL to get reflectance from 4 fluxes
  Ref <- PRO4SAIL(Spec_Sensor = SpecPROSPECT_Sensor, Input_PROSPECT = InVar,
                  TypeLidf = TypeLidf, LIDFa = InVar$LIDFa, LIDFb = InVar$LIDFb,
                  lai = InVar$lai, q = InVar$q, tts = InVar$tts,
                  tto = InVar$tto, psi = InVar$psi, rsoil = rsoil)
  # Computes bidirectional reflectance factor based on outputs from PROSAIL and sun position
  brfMOD <- Compute_BRF(rdot = Ref$rdot, rsot = Ref$rsot,
                        tts = InVar$tts, SpecATM_Sensor = SpecATM_Sensor)
  # compute cost
  fc <- CostVal_RMSE_PROSAIL(brfMES = brfMES, brfMOD$BRF, xprior,
                             PriorInfoMean = PriorInfoMean,
                             PriorInfoSD = PriorInfoSD,
                             WeightPrior = WeightPrior)
  return(fc)
}


#' Value of the cost criterion to minimize during PROSAIL inversion
#' @param brfMES numeric. Measured bidirectional reflectance
#' @param brfMOD numeric. Simulated bidirectional reflectance
#' @param xprior list. values of the parameters for which prior information is provided
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard deviation of parameters defined as xprior
#' @param WeightPrior numeric. Weight to be applied on prior information to modulate its importance
#'
#' @return res list. Includes Parms2Estimate, Parms2Set,
#' InitialGuess, LowerBound, UpperBound, ParmSet, InVar
#' @export
CostVal_RMSE_PROSAIL  <- function(brfMES, brfMOD, xprior, PriorInfoMean = NULL,
                                  PriorInfoSD = NULL, WeightPrior = 0.01){

  fc <- sqrt(sum((brfMES-brfMOD)**2)/length(brfMES))
  if (!is.null(PriorInfoMean)){
    fc <- fc + WeightPrior*mean(as.numeric((xprior-PriorInfoMean)/PriorInfoSD)**2)
  }
  return(fc)
}

#' function identifying which parameters should be estimated during inversion
#'
#' @param InitialGuess list. Initial values during inversion
#' @param LowerBound list. Lower bound values during inversion
#' @param UpperBound list. Upper bound values during inversion
#' @param ParmSet list. Values for variables out of inversion
#'
#' @return res list. includes
#' - Parms2Estimate
#' - Parms2Set
#' - InitialGuess
#' - LowerBound
#' - UpperBound
#' - ParmSet
#' - InVar
#' fc estimates of the parameters
#' @export
WhichParameters2Invert <- function(InitialGuess, LowerBound,
                                   UpperBound, ParmSet) {

  # define all parameters which can be assessed through iterativbe optimization
  InVar <- data.frame('CHL' = 0,'CAR' = 0,'ANT' = 0,'BROWN' = 0,'EWT' = 0,
                      'LMA' = 0,'PROT' = 0,'CBC' = 0,'N' = 0,'alpha' = 40,
                      'LIDFa' = 0,'LIDFb' = 0,'lai' = 0,'q' = 0,
                      'tts' = 0,'tto' = 0,'psi' = 0,'psoil' = 0)
  AllParms <- names(InVar)
  Parms2Estimate <- ParmSet_Final <- c()
  InitialGuess_Update <- LowerBound_Update <- UpperBound_Update <-
    ParmSet_Update <- list()

  # set parameters to user value defined in ParmSet
  for (parm in AllParms){
    if (parm %in% names(InitialGuess) & !parm %in% names(ParmSet)){
      Parms2Estimate <- c(Parms2Estimate,parm)
      InitialGuess_Update[[parm]] <- InitialGuess[[parm]]
      LowerBound_Update[[parm]] <- LowerBound[[parm]]
      UpperBound_Update[[parm]] <- UpperBound[[parm]]
    } else if (!parm %in% names(InitialGuess) & parm %in% names(ParmSet)){
      ParmSet_Final <- c(ParmSet_Final,parm)
      ParmSet_Update[[parm]] <- ParmSet[[parm]]
    }
  }
  InitialGuess_Update <- data.frame(InitialGuess_Update)
  LowerBound_Update <- data.frame(LowerBound_Update)
  UpperBound_Update <- data.frame(UpperBound_Update)
  ParmSet_Update <- data.frame(ParmSet_Update)
  return(list('Parms2Estimate' = Parms2Estimate,
              'Parms2Set' = ParmSet_Final,
              'InitialGuess' = InitialGuess_Update,
              'LowerBound' = LowerBound_Update,
              'UpperBound' = UpperBound_Update,
              'ParmSet' = ParmSet_Update,
              'InVar' = InVar))
}

#' function identifying for which parameters prior information is known
#'
#' @param PriorInfoMean list. mean value expected for a list of parameters
#' @param PriorInfoSD list. dtandard deviation to the mean expected for a list of parameters
#'
#' @return fc estimates of the parameters
#' @export
WhichParmPrior <- function(PriorInfoMean, PriorInfoSD) {

  Parms2Prior <- c()
  listParms <- c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'PROT', 'CBC', 'N',
                 'alpha', 'LIDFa', 'LIDFb', 'lai', 'q', 'tts', 'tto', 'psi',
                 'psoil')
  PriorInfoMean_Update <- PriorInfoSD_Update <- list()
  for (parm in listParms){
    if (parm %in% names(PriorInfoMean) & parm %in% names(PriorInfoSD)){
      Parms2Prior <- c(Parms2Prior,parm)
      PriorInfoMean_Update[[parm]] <- PriorInfoMean[[parm]]
      PriorInfoSD_Update[[parm]] <- PriorInfoSD[[parm]]
    }
  }
  PriorInfoMean_Update <- data.frame(PriorInfoMean_Update)
  PriorInfoSD_Update <- data.frame(PriorInfoSD_Update)
  return(list('Parms2Prior' = Parms2Prior,
              'PriorInfoMean' = PriorInfoMean_Update,
              'PriorInfoSD' = PriorInfoSD_Update))
}
