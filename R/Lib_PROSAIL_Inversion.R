# ==============================================================================
# prosail
# Lib_PROSAIL_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to PROSAIL inversion using iterative
# optimization
# ==============================================================================

#' Performs PROSAIL inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param brfMES numeric. measured BRF to be used during inversion
#' @param InitialGuess list. Initial guess of Parms2Estimate for the inversion
#' @param LowerBound list. Lower bound
#' @param UpperBound list. Upper bound
#' @param SpecPROSPECT_Sensor  list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
#' @param SpecSOIL_Sensor list. includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param TypeLidf numeric. Type of leaf inclination distribution function
#' @param ParmSet  list. Parameters set to a fixed value by user
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard deviation of parameters defined as xprior
#' @param WeightPrior numeric. Weight to be applied on prior information to modulate its importance
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom pracma fmincon
#' @export

Invert_PROSAIL  <- function(brfMES,InitialGuess,LowerBound,UpperBound,
                            SpecPROSPECT_Sensor,SpecATM_Sensor,SpecSOIL_Sensor,TypeLidf,
                            ParmSet,MeritFunction = 'Merit_RMSE_PROSAIL',PriorInfoMean=NULL,PriorInfoSD=NULL,WeightPrior = 0.01){

  # define parameters included in inversion or set
  ParmInv <- WhichParameters2Invert(InitialGuess,LowerBound,UpperBound,ParmSet)
  if (!is.null(PriorInfoMean) & !is.null(PriorInfoSD)){
    PriorInfo <- WhichParmPRior(PriorInfoMean,PriorInfoSD)
    PriorInfoMean <- PriorInfo$PriorInfoMean
    PriorInfoSD <- PriorInfo$PriorInfoSD
    Parms2Prior <- PriorInfo$Parms2Prior
  }
  ParmInv$InVar[ParmInv$Parms2Set] <- ParmInv$ParmSet

  # update init value and lower/upper boundaries for inversion based on Vars2Estimate
  xinit = as.vector(ParmInv$InitialGuess,mode='numeric')
  lb    = as.vector(ParmInv$LowerBound,mode='numeric')
  ub    = as.vector(ParmInv$UpperBound,mode='numeric')
  resInv   = fmincon(x0 = xinit, fn = MeritFunction, gr = NULL,brfMES = brfMES,
                  SpecPROSPECT_Sensor=SpecPROSPECT_Sensor, SpecSOIL_Sensor=SpecSOIL_Sensor, SpecATM_Sensor=SpecATM_Sensor,
                  Parms2Estimate=ParmInv$Parms2Estimate, Parm2Set=ParmInv$Parms2Set, ParmSet=ParmInv$ParmSet ,
                  InVar=ParmInv$InVar ,TypeLidf=TypeLidf,
                  PriorInfoMean =PriorInfoMean,PriorInfoSD =PriorInfoSD,Parms2Prior=Parms2Prior,WeightPrior=WeightPrior,
                  method = "SQP",A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                  lb = lb, ub = ub, hin = NULL, heq = NULL,tol = 1e-08,
                  maxfeval = 2000, maxiter = 2000)

  ParmInv$InVar[ParmInv$Parms2Estimate] <- resInv$par
  return(ParmInv$InVar)
}

#' Merit function for PROSAIL inversion
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param brfMES  numeric. measured BRF
#' @param SpecPROSPECT_Sensor list. Includes optical constants for PROSPECT
#' @param SpecSOIL_Sensor list. Includes soil reflectance (either 2 references for optimization of psoil, or a unique spectrun)
#' @param SpecATM_Sensor list. Includes diffuse and direct light
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Parms2Estimate  numeric. rank of variables to be inverted
#' to be estimated through inversion
#' @param Parm2Set  numeric. rank of variables to be set out of inversion
#' @param ParmSet  list. value of variables to be set out of inversion
#' @param InVar dataframe. full set of PROSAIL input variables
#' @param TypeLidf  numeric. tyoe of leaf inclination distribution function
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard deviation of parameters defined as xprior
#' @param Parms2Prior numeric. rank of the parameters which should be used with prior information
#' @param WeightPrior numeric. Weight to be applied on prior information to modulate its importance
#'
#' @return fc estimates of the parameters
#' @export
Merit_RMSE_PROSAIL <- function(xinit,brfMES,SpecPROSPECT_Sensor,SpecSOIL_Sensor,SpecATM_Sensor,
                               Parms2Estimate,Parm2Set,ParmSet,InVar,TypeLidf,
                               PriorInfoMean=NULL,PriorInfoSD=NULL,Parms2Prior=NULL,WeightPrior=0.01){

  xinit[xinit<0] = 0
  InVar[Parms2Estimate] <- xinit
  xprior <-InVar[Parms2Prior]
  rsoil <- InVar$psoil*SpecSOIL_Sensor$Dry_Soil+(1-InVar$psoil)*SpecSOIL_Sensor$Wet_Soil
  # call PROSAIL to get reflectance from 4 fluxes
  Ref <- PRO4SAIL(SpecPROSPECT_Sensor,CHL = InVar$CHL, CAR = InVar$CAR, ANT = InVar$ANT,
                  EWT = InVar$EWT, LMA = InVar$LMA, N = InVar$N,
                  TypeLidf = TypeLidf,LIDFa = InVar$LIDFa,LIDFb = InVar$LIDFb,lai = InVar$lai,
                  q = InVar$q,tts = InVar$tts,tto = InVar$tto,psi = InVar$psi,rsoil = rsoil)
  # Computes bidirectional reflectance factor based on outputs from PROSAIL and sun position
  brfMOD <- Compute_BRF(Ref$rdot,Ref$rsot,InVar$tts,SpecATM_Sensor)
  # compute cost
  fc <- CostVal_RMSE_PROSAIL(brfMES,brfMOD,xprior,PriorInfoMean=PriorInfoMean,PriorInfoSD=PriorInfoSD,WeightPrior=WeightPrior)
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
CostVal_RMSE_PROSAIL  <- function(brfMES,brfMOD,xprior,PriorInfoMean=NULL,PriorInfoSD=NULL,WeightPrior=0.01){

  fc = sqrt(sum((brfMES-brfMOD)**2)/length(brfMES))
  if (!is.null(PriorInfoMean)){
    # message('fc')
    # print(fc)
    fc = fc + WeightPrior*mean(as.numeric((xprior-PriorInfoMean)/PriorInfoSD)**2)
    # message('prior')
    # print(0.01*mean(as.numeric((xprior-PriorInfoMean)/PriorInfoSD)**2))
    # print(PriorInfoMean)
    # print(PriorInfoSD)
    # print((((xprior-PriorInfoMean)/PriorInfoSD)**2))
    # print(mean(((xprior-PriorInfoMean)/PriorInfoSD)**2))
    # print(0.01*mean(((xprior-PriorInfoMean)/PriorInfoSD)**2))
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
WhichParameters2Invert <- function(InitialGuess,LowerBound,UpperBound,ParmSet) {

  Parms2Estimate <- c()
  ParmSet_Final <- c()
  InitialGuess_Update <- LowerBound_Update <- UpperBound_Update <- ParmSet_Update <- data.frame(row.names = c('Val'))

  # set parameters to user value defined in ParmSet
  if ('CHL'%in%names(InitialGuess) & !'CHL'%in%names(ParmSet)){
    Parms2Estimate <- c(Parms2Estimate,1)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'CHL' = InitialGuess$CHL)
    LowerBound_Update <- data.frame(LowerBound_Update,'CHL' = LowerBound$CHL)
    UpperBound_Update <- data.frame(UpperBound_Update,'CHL' = UpperBound$CHL)
  } else if (!'CHL'%in%names(InitialGuess) & 'CHL'%in%names(ParmSet)){
    ParmSet_Final <- c(ParmSet_Final,1)
    ParmSet_Update <- data.frame(ParmSet_Update,'CHL' = ParmSet$CHL)
  }
  if ('CAR'%in%names(InitialGuess) & !'CAR'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,2)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'CAR' = InitialGuess$CAR)
    LowerBound_Update <- data.frame(LowerBound_Update,'CAR' = LowerBound$CAR)
    UpperBound_Update <- data.frame(UpperBound_Update,'CAR' = UpperBound$CAR)
  } else if (!'CAR'%in%names(InitialGuess) & 'CAR'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,2)
    ParmSet_Update <- data.frame(ParmSet_Update,'CAR' = ParmSet$CAR)
  }
  if ('ANT'%in%names(InitialGuess) & !'ANT'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,3)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'ANT' = InitialGuess$ANT)
    LowerBound_Update <- data.frame(LowerBound_Update,'ANT' = LowerBound$ANT)
    UpperBound_Update <- data.frame(UpperBound_Update,'ANT' = UpperBound$ANT)
  } else if (!'ANT'%in%names(InitialGuess) & 'ANT'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,3)
    ParmSet_Update <- data.frame(ParmSet_Update,'ANT' = ParmSet$ANT)
  }
  if ('BROWN'%in%names(InitialGuess) & !'BROWN'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,4)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'BROWN' = InitialGuess$BROWN)
    LowerBound_Update <- data.frame(LowerBound_Update,'BROWN' = LowerBound$BROWN)
    UpperBound_Update <- data.frame(UpperBound_Update,'BROWN' = UpperBound$BROWN)
  } else if (!'BROWN'%in%names(InitialGuess) & 'BROWN'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,4)
    ParmSet_Update <- data.frame(ParmSet_Update,'BROWN' = ParmSet$BROWN)
  }
  if ('EWT'%in%names(InitialGuess) & !'EWT'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,5)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'EWT' = InitialGuess$EWT)
    LowerBound_Update <- data.frame(LowerBound_Update,'EWT' = LowerBound$EWT)
    UpperBound_Update <- data.frame(UpperBound_Update,'EWT' = UpperBound$EWT)
  } else if (!'EWT'%in%names(InitialGuess) & 'EWT'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,5)
    ParmSet_Update <- data.frame(ParmSet_Update,'EWT' = ParmSet$EWT)
  }
  if ('LMA'%in%names(InitialGuess) & !'LMA'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,6)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'LMA' = InitialGuess$LMA)
    LowerBound_Update <- data.frame(LowerBound_Update,'LMA' = LowerBound$LMA)
    UpperBound_Update <- data.frame(UpperBound_Update,'LMA' = UpperBound$LMA)
  } else if (!'LMA'%in%names(InitialGuess) & 'LMA'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,6)
    ParmSet_Update <- data.frame(ParmSet_Update,'LMA' = ParmSet$LMA)
  }
  if ('PROT'%in%names(InitialGuess) & !'PROT'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,7)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'PROT' = InitialGuess$PROT)
    LowerBound_Update <- data.frame(LowerBound_Update,'PROT' = LowerBound$PROT)
    UpperBound_Update <- data.frame(UpperBound_Update,'PROT' = UpperBound$PROT)
  } else if (!'PROT'%in%names(InitialGuess) & 'PROT'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,7)
    ParmSet_Update <- data.frame(ParmSet_Update,'PROT' = ParmSet$PROT)
  }
  if ('CBC'%in%names(InitialGuess) & !'CBC'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,8)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'CBC' = InitialGuess$CBC)
    LowerBound_Update <- data.frame(LowerBound_Update,'CBC' = LowerBound$CBC)
    UpperBound_Update <- data.frame(UpperBound_Update,'CBC' = UpperBound$CBC)
  } else if (!'CBC'%in%names(InitialGuess) & 'CBC'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,8)
    ParmSet_Update <- data.frame(ParmSet_Update,'CBC' = ParmSet$CBC)
  }
  if ('N'%in%names(InitialGuess) & !'N'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,9)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'N' = InitialGuess$N)
    LowerBound_Update <- data.frame(LowerBound_Update,'N' = LowerBound$N)
    UpperBound_Update <- data.frame(UpperBound_Update,'N' = UpperBound$N)
  } else if (!'N'%in%names(InitialGuess) & 'N'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,9)
    ParmSet_Update <- data.frame(ParmSet_Update,'N' = ParmSet$N)
  }
  if ('alpha'%in%names(InitialGuess) & !'alpha'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,10)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'alpha' = InitialGuess$alpha)
    LowerBound_Update <- data.frame(LowerBound_Update,'alpha' = LowerBound$alpha)
    UpperBound_Update <- data.frame(UpperBound_Update,'alpha' = UpperBound$alpha)
  } else if (!'alpha'%in%names(InitialGuess) & 'alpha'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,10)
    ParmSet_Update <- data.frame(ParmSet_Update,'alpha' = ParmSet$alpha)
  }
  if ('LIDFa'%in%names(InitialGuess) & !'LIDFa'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,11)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'LIDFa' = InitialGuess$LIDFa)
    LowerBound_Update <- data.frame(LowerBound_Update,'LIDFa' = LowerBound$LIDFa)
    UpperBound_Update <- data.frame(UpperBound_Update,'LIDFa' = UpperBound$LIDFa)
  } else if (!'LIDFa'%in%names(InitialGuess) & 'LIDFa'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,11)
    ParmSet_Update <- data.frame(ParmSet_Update,'LIDFa' = ParmSet$LIDFa)
  }
  if ('LIDFb'%in%names(InitialGuess) & !'LIDFb'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,12)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'LIDFb' = InitialGuess$LIDFb)
    LowerBound_Update <- data.frame(LowerBound_Update,'LIDFb' = LowerBound$LIDFb)
    UpperBound_Update <- data.frame(UpperBound_Update,'LIDFb' = UpperBound$LIDFb)
  } else if (!'LIDFb'%in%names(InitialGuess) & 'LIDFb'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,12)
    ParmSet_Update <- data.frame(ParmSet_Update,'LIDFb' = ParmSet$LIDFb)
  }
  if ('lai'%in%names(InitialGuess) & !'lai'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,13)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'lai' = InitialGuess$lai)
    LowerBound_Update <- data.frame(LowerBound_Update,'lai' = LowerBound$lai)
    UpperBound_Update <- data.frame(UpperBound_Update,'lai' = UpperBound$lai)
  } else if (!'lai'%in%names(InitialGuess) & 'lai'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,13)
    ParmSet_Update <- data.frame(ParmSet_Update,'lai' = ParmSet$lai)
  }
  if ('q'%in%names(InitialGuess) & !'q'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,14)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'q' = InitialGuess$q)
    LowerBound_Update <- data.frame(LowerBound_Update,'q' = LowerBound$q)
    UpperBound_Update <- data.frame(UpperBound_Update,'q' = UpperBound$q)
  } else if (!'q'%in%names(InitialGuess) & 'q'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,14)
    ParmSet_Update <- data.frame(ParmSet_Update,'q' = ParmSet$q)
  }
  if ('tts'%in%names(InitialGuess) & !'tts'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,15)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'tts' = InitialGuess$tts)
    LowerBound_Update <- data.frame(LowerBound_Update,'tts' = LowerBound$tts)
    UpperBound_Update <- data.frame(UpperBound_Update,'tts' = UpperBound$tts)
  } else if (!'tts'%in%names(InitialGuess) & 'tts'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,15)
    ParmSet_Update <- data.frame(ParmSet_Update,'tts' = ParmSet$tts)
  }
  if ('tto'%in%names(InitialGuess) & !'tto'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,16)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'tto' = InitialGuess$tto)
    LowerBound_Update <- data.frame(LowerBound_Update,'tto' = LowerBound$tto)
    UpperBound_Update <- data.frame(UpperBound_Update,'tto' = UpperBound$tto)
  } else if (!'tto'%in%names(InitialGuess) & 'tto'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,16)
    ParmSet_Update <- data.frame(ParmSet_Update,'tto' = ParmSet$tto)
  }
  if ('psi'%in%names(InitialGuess) & !'psi'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,17)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'psi' = InitialGuess$psi)
    LowerBound_Update <- data.frame(LowerBound_Update,'psi' = LowerBound$psi)
    UpperBound_Update <- data.frame(UpperBound_Update,'psi' = UpperBound$psi)
  } else if (!'psi'%in%names(InitialGuess) & 'psi'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,17)
    ParmSet_Update <- data.frame(ParmSet_Update,'psi' = ParmSet$psi)
  }
  if ('psoil'%in%names(InitialGuess) & !'psoil'%in%names(ParmSet)){
    Parms2Estimate = c(Parms2Estimate,18)
    InitialGuess_Update <- data.frame(InitialGuess_Update,'psoil' = InitialGuess$psoil)
    LowerBound_Update <- data.frame(LowerBound_Update,'psoil' = LowerBound$psoil)
    UpperBound_Update <- data.frame(UpperBound_Update,'psoil' = UpperBound$psoil)
  } else if (!'psoil'%in%names(InitialGuess) & 'psoil'%in%names(ParmSet)){
    ParmSet_Final = c(ParmSet_Final,18)
    ParmSet_Update <- data.frame(ParmSet_Update,'psoil' = ParmSet$psoil)
  }
  InVar <- data.frame('CHL'=0,'CAR'=0,'ANT'=0,'BROWN'=0,'EWT'=0,
                      'LMA'=0,'PROT'=0,'CBC'=0,'N'=0,'alpha'=0,
                      'LIDFa'=0,'LIDFb'=0,'lai'=0,'q'=0,
                      'tts'=0,'tto'=0,'psi'=0,'psoil'=0)

  res <- list('Parms2Estimate'=Parms2Estimate,'Parms2Set'=ParmSet_Final,
             'InitialGuess'=InitialGuess_Update,'LowerBound'= LowerBound_Update,
             'UpperBound'=UpperBound_Update,'ParmSet'= ParmSet_Update,'InVar'= InVar)
  return(res)
}

#' function identifying for which parameters prior information is known
#'
#' @param PriorInfoMean list. mean value expected for a list of parameters
#' @param PriorInfoSD list. dtandard deviation to the mean expected for a list of parameters
#'
#' @return fc estimates of the parameters
#' @export
WhichParmPRior <- function(PriorInfoMean,PriorInfoSD) {

  Parms2Prior <- c()
  PriorInfoMean_Update <- data.frame(row.names = c('Val'))
  PriorInfoSD_Update <- data.frame(row.names = c('Val'))

  # set parameters to user value defined in ParmSet
  if ('CHL'%in%names(PriorInfoMean) & 'CHL'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,1)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'CHL' = PriorInfoMean$CHL)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'CHL' = PriorInfoSD$CHL)
  }
  if ('CAR'%in%names(PriorInfoMean) & 'CAR'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,2)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'CAR' = PriorInfoMean$CAR)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'CAR' = PriorInfoSD$CAR)
  }
  if ('ANT'%in%names(PriorInfoMean) & 'ANT'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,3)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'ANT' = PriorInfoMean$ANT)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'ANT' = PriorInfoSD$ANT)
  }
  if ('BROWN'%in%names(PriorInfoMean) & 'BROWN'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,4)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'BROWN' = PriorInfoMean$BROWN)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'BROWN' = PriorInfoSD$BROWN)
  }
  if ('EWT'%in%names(PriorInfoMean) & 'EWT'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,5)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'EWT' = PriorInfoMean$EWT)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'EWT' = PriorInfoSD$EWT)
  }
  if ('LMA'%in%names(PriorInfoMean) & 'LMA'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,6)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'LMA' = PriorInfoMean$LMA)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'LMA' = PriorInfoSD$LMA)
  }
  if ('PROT'%in%names(PriorInfoMean) & 'PROT'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,7)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'PROT' = PriorInfoMean$PROT)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'PROT' = PriorInfoSD$PROT)
  }
  if ('CBC'%in%names(PriorInfoMean) & 'CBC'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,8)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'CBC' = PriorInfoMean$CBC)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'CBC' = PriorInfoSD$CBC)
  }
  if ('N'%in%names(PriorInfoMean) & 'N'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,9)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'N' = PriorInfoMean$N)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'N' = PriorInfoSD$N)
  }
  if ('alpha'%in%names(PriorInfoMean) & 'alpha'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,10)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'alpha' = PriorInfoMean$alpha)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'alpha' = PriorInfoSD$alpha)
  }
  if ('LIDFa'%in%names(PriorInfoMean) & 'LIDFa'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,11)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'LIDFa' = PriorInfoMean$LIDFa)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'LIDFa' = PriorInfoSD$LIDFa)
  }
  if ('LIDFb'%in%names(PriorInfoMean) & 'LIDFb'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,12)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'LIDFb' = PriorInfoMean$LIDFb)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'LIDFb' = PriorInfoSD$LIDFb)
  }
  if ('lai'%in%names(PriorInfoMean) & 'lai'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,13)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'lai' = PriorInfoMean$lai)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'lai' = PriorInfoSD$lai)
  }
  if ('q'%in%names(PriorInfoMean) & 'q'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,14)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'q' = PriorInfoMean$q)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'q' = PriorInfoSD$q)
  }
  if ('tts'%in%names(PriorInfoMean) & 'tts'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,15)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'tts' = PriorInfoMean$tts)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'tts' = PriorInfoSD$tts)
  }
  if ('tto'%in%names(PriorInfoMean) & 'tto'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,16)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'tto' = PriorInfoMean$tto)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'tto' = PriorInfoSD$tto)
  }
  if ('psi'%in%names(PriorInfoMean) & 'psi'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,17)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'psi' = PriorInfoMean$psi)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'psi' = PriorInfoSD$psi)
  }
  if ('psoil'%in%names(PriorInfoMean) & 'psoil'%in%names(PriorInfoSD)){
    Parms2Prior <- c(Parms2Prior,18)
    PriorInfoMean_Update <- data.frame(PriorInfoMean_Update,'psoil' = PriorInfoMean$psoil)
    PriorInfoSD_Update <- data.frame(PriorInfoSD_Update,'psoil' = PriorInfoSD$psoil)
  }
  res = list('Parms2Prior'=Parms2Prior,'PriorInfoMean'=PriorInfoMean_Update,'PriorInfoSD'=PriorInfoSD_Update)
  return(res)
}
