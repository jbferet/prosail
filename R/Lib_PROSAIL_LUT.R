# ============================================================================= =
# prosail
# Lib_PROSAIL_LUT.R
# ============================================================================= =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================= =
# This Library includes functions dedicated to generating PROSAIL LUTs
# ============================================================================= =

#' This function provides a LUT corresponding to PROSAIL input parameters to
#' reproduce the distribution of parameters used in the ATBD
#'
#' @param nbSamples numeric. number of samples to generate
#' @param GeomAcq list. min and max values for tts, tto and psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#'
#' @return InputPROSAIL list. list of PROSAIL input parameters
#' @importFrom truncnorm rtruncnorm
#' @export

get_atbd_LUT_input <- function(nbSamples = 2000, GeomAcq = NULL, Codist_LAI = TRUE){

  # define paremertization for truncated gaussians
  TgaussParms <- list()
  TgaussParms$min <- data.frame('lai' = 0, 'LIDFa' = 30, 'q' = 0.1, 'N' = 1.2, 'CHL' = 20, 'LMA' = 0.003, 'Cw_rel' = 0.6, 'BROWN' = 0.0, 'psoil' = 0)
  TgaussParms$max <- data.frame('lai' = 15, 'LIDFa' = 80, 'q' = 0.5, 'N' = 1.8, 'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.85, 'BROWN' = 2.0, 'psoil' = 1)
  TgaussParms$mean <- data.frame('lai' = 2, 'LIDFa' = 60, 'q' = 0.2, 'N' = 1.5, 'CHL' = 45, 'LMA' = 0.005, 'Cw_rel' = 0.75, 'BROWN' = 0.0, 'psoil' = 0.25)
  TgaussParms$sd <- data.frame('lai' = 3, 'LIDFa' = 30, 'q' = 0.5, 'N' = 0.3, 'CHL' = 30, 'LMA' = 0.005, 'Cw_rel' = 0.08, 'BROWN' = 0.30, 'psoil' = 0.6)

  # get distribution corresponding to gaussians
  InputPROSAIL <- list()
  for (parm in names(TgaussParms$min)){
    set.seed(42)
    InputPROSAIL[[parm]] <- truncnorm::rtruncnorm(n = nbSamples,
                                                  a = TgaussParms$min[[parm]], b = TgaussParms$max[[parm]],
                                                  mean = TgaussParms$mean[[parm]], sd = TgaussParms$sd[[parm]])
  }
  InputPROSAIL <- data.frame(InputPROSAIL)

  # define co-distribution with LAI
  if (Codist_LAI==TRUE){
    Codist_LAI <- list()
    Codist_LAI$Vmin0 <- data.frame('LIDFa' = 30, 'q' = 0.1, 'N' = 1.2, 'CHL' = 20, 'LMA' = 0.003, 'Cw_rel' = 0.6, 'BROWN' = 0.0, 'psoil' = 0)
    Codist_LAI$Vmax0 <- data.frame('LIDFa' = 80, 'q' = 0.5, 'N' = 1.8, 'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.85, 'BROWN' = 2.0, 'psoil' = 1)
    Codist_LAI$VminLAImax <- data.frame('LIDFa' = 55, 'q' = 0.1, 'N' = 1.3, 'CHL' = 45, 'LMA' = 0.005, 'Cw_rel' = 0.70, 'BROWN' = 0.0, 'psoil' = 0)
    Codist_LAI$VmaxLAImax <- data.frame('LIDFa' = 65, 'q' = 0.5, 'N' = 1.8, 'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.80, 'BROWN' = 0.2, 'psoil' = 0.4)
    for (parm in names(Codist_LAI$Vmin0)){
      Vstar <- get_codistributions(V = InputPROSAIL[[parm]],
                                   LAI = InputPROSAIL$lai, MaxLAI = TgaussParms$max$lai,
                                   Vmin0 = Codist_LAI$Vmin0[[parm]],
                                   Vmax0 = Codist_LAI$Vmax0[[parm]],
                                   VminLAImax = Codist_LAI$VminLAImax[[parm]],
                                   VmaxLAImax = Codist_LAI$VmaxLAImax[[parm]])
      InputPROSAIL[[parm]] <- Vstar
    }
  }

  # convert CW_rel into EWT
  InputPROSAIL$EWT <- ((InputPROSAIL$LMA)/(1-InputPROSAIL$Cw_rel))-InputPROSAIL$LMA
  # set ANT, PROT and CBC to 0
  InputPROSAIL$ANT <- InputPROSAIL$PROT <- InputPROSAIL$CBC <- 0
  # set CAR to 0.25*CHL
  InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL
  # set geometry of acquisition
  if (is.null(GeomAcq)){
    GeomAcq <- list()
    GeomAcq$min <- GeomAcq$max <- list()
    GeomAcq$min$tto <- 0
    GeomAcq$max$tto <- 10
    GeomAcq$min$tts <- 20
    GeomAcq$max$tts <- 30
    GeomAcq$min$psi <- 0
    GeomAcq$max$psi <- 360

  }
  InputPROSAIL$tts <- runif(n = nbSamples, min = GeomAcq$min$tts,max=GeomAcq$max$tts)
  InputPROSAIL$tto <- runif(n = nbSamples, min = GeomAcq$min$tto,max=GeomAcq$max$tto)
  InputPROSAIL$psi <- runif(n = nbSamples, min = GeomAcq$min$psi,max=GeomAcq$max$psi)
  return(InputPROSAIL)
}

#' This function adjusts variable values based on co-distribution rules as defined in ATBD
#' co-distributions are all related to LAI
#'
#' @param V numeric.
#' @param LAI numeric.
#' @param MaxLAI numeric.
#' @param Vmin0 numeric.
#' @param Vmax0 numeric.
#' @param VminLAImax numeric.
#' @param VmaxLAImax numeric.
#'
#' @return Vstar numeric.
#' @export

get_codistributions <- function(V, LAI, MaxLAI, Vmin0, Vmax0, VminLAImax, VmaxLAImax){

  VminLAI <- Vmin0 + (LAI*(VminLAImax-Vmin0)/MaxLAI)
  VmaxLAI <- Vmax0 + (LAI*(VmaxLAImax-Vmax0)/MaxLAI)
  Vstar <- VminLAI+((V-Vmin0)*(VmaxLAI-VminLAI)/(Vmax0-Vmin0))
  return(Vstar)
}


#' This function sets default values for PROSAIL LUT simulation when not defined by user
#'
#' @param TypeDistrib list. specify if uniform or Gaussian distribution to be applied. default = Uniform
#' @param GaussianDistrib  list. Mean value and STD corresponding to the parameters sampled with gaussian distribution
#' @param minval list. Defines the minimum value to be set for a list of parameters randomly produced
#' @param maxval list. Defines the maximum value to be set for a list of parameters randomly produced
#'
#' @return res list. list of default values corresponding to NULL input parameters
#' @export

get_default_LUT_input <- function(TypeDistrib = NULL,
                                  GaussianDistrib = NULL,
                                  minval = NULL,
                                  maxval = NULL){
  # define uniform / gaussian distribution
  if (is.null(TypeDistrib)) TypeDistrib <- data.frame('CHL'='Uniform', 'CAR'='Uniform', 'ANT' = 'Uniform', 'BROWN'='Uniform',
                                                      'EWT' = 'Uniform', 'LMA' = 'Uniform', 'N' = 'Uniform',
                                                      'psoil' = 'Uniform', 'LIDFa' = 'Uniform', 'lai' = 'Uniform', 'q'='Uniform',
                                                      'tto' = 'Uniform','tts' = 'Uniform', 'psi' = 'Uniform')
  # define mean / sd for gaussian
  if (is.null(GaussianDistrib)) GaussianDistrib <- list('Mean'=NULL,'Std'=NULL)
  # define min and max values
  if (is.null(minval)) minval <- data.frame('CHL'=10,'CAR'=0,'EWT' = 0.01,'ANT' = 0,'LMA' = 0.005,'N' = 1.0,'psoil' = 0.0, 'BROWN'=0.0,
                                            'LIDFa' = 20, 'lai' = 0.5,'q'=0.1,'tto' = 0,'tts' = 20, 'psi' = 80)
  if (is.null(maxval)) maxval <- data.frame('CHL'=75,'CAR'=15,'EWT' = 0.03,'ANT' = 2,'LMA' = 0.03,'N' = 2.0, 'psoil' = 1.0, 'BROWN'=0.5,
                                            'LIDFa' = 70, 'lai' = 7,'q'=0.2,'tto' = 5,'tts' = 30, 'psi' = 110)
  res <- list('TypeDistrib' = TypeDistrib,
              'GaussianDistrib' = GaussianDistrib,
              'minval' = minval, 'maxval' = maxval)
  return(res)
}

#' This function generates distribution of biophysical parameters used as input parameters in PRO4SAIL
#'
#' @param minval list. Defines the minimum value to be set for a list of parameters randomly produced
#' @param maxval list. Defines the maximum value to be set for a list of parameters randomly produced
#' @param ParmSet list.Defines the parameters to be set to a given value
#' @param nbSamples numeric. Number of samples to be generated
#' @param TypeDistrib list. specify if uniform or Gaussian distribution to be applied. default = Uniform
#' @param Mean list. mean value for parameters with Gaussian distribution
#' @param Std list. standard deviation for parameters with Gaussian distribution
#' @param Force4LowLAI boolean. set to TRUE if needed to force distribution of leaf chemistry to low values for low LAI
#'
#' @return InputPROSAIL list. list of nbSamples input parameters for PRO4SAIL simulation
#' @importFrom stats runif rnorm sd
#' @export
get_distribution_input_prosail <- function(minval,maxval,ParmSet,nbSamples,
                                           TypeDistrib = data.frame('CHL'='Uniform',
                                                                    'CAR'='Uniform',
                                                                    'EWT' = 'Uniform',
                                                                    'ANT' = 'Uniform',
                                                                    'LMA' = 'Uniform',
                                                                    'PROT' = 'Uniform',
                                                                    'CBC' = 'Uniform',
                                                                    'N' = 'Uniform',
                                                                    'BROWN'='Uniform',
                                                                    'psoil' = 'Uniform',
                                                                    'LIDFa' = 'Uniform',
                                                                    'LIDFb' = 'Uniform',
                                                                    'lai' = 'Uniform',
                                                                    'q'='Uniform',
                                                                    'tto' = 'Uniform',
                                                                    'tts' = 'Uniform',
                                                                    'psi' = 'Uniform',
                                                                    'alpha' = 'Uniform'),
                                           Mean = NULL, Std = NULL, Force4LowLAI=TRUE){

  # define all input parameters from PROSAIL
  InVar <- c('CHL','CAR','ANT','BROWN','EWT','LMA',
             'PROT','CBC','N','alpha','LIDFa','LIDFb',
             'lai','q','tts','tto','psi','psoil','TypeLidf')

  InputPROSAIL <- list('CHL'=c(),'CAR'=c(),'ANT'=c(),'BROWN'=c(),'EWT'=c(),
                       'LMA'=c(),'PROT'=c(),'CBC'=c(),'N'=c(),'alpha'=c(),
                       'LIDFa'=c(),'LIDFb'=c(),'lai'=c(),'q'=c(),
                       'tts'=c(),'tto'=c(),'psi'=c(),'psoil'=c(),'TypeLidf'=c())
  Default <- data.frame('CHL'=0,'CAR'=0,'ANT'=0,'BROWN'=0,'EWT'=0,
                        'LMA'=0,'PROT'=0,'CBC'=0,'N'=1.5,'alpha'=40,
                        'LIDFa'=0,'LIDFb'=0,'lai'=2,'q'=0,
                        'tts'=0,'tto'=0,'psi'=0,'psoil'=1,'TypeLidf'=2)

  # which input parameters should be randomly sampled?
  ParmRand <- which(is.element(InVar,names(minval))==TRUE)
  # which input parameters should be set to a constant value?
  ParmCte <- which(is.element(InVar,names(ParmSet))==TRUE)
  # if some parameters are defined to be sampled randomly and set to constant value
  if (length(intersect(ParmRand,ParmCte))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(InVar[intersect(ParmRand,ParmCte)])
    message('will be randomly sampled based on the min and max values defined by minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  AllParm <- c(ParmRand,ParmCte)
  FullList <- seq(1,length(InVar),by=1)
  Set2Default <- c()
  if (length(setdiff(FullList,AllParm))>0){
    message('WARNING: some elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(InVar[setdiff(FullList,AllParm)])
    message('These parameters will be set to their default value')
    print(Default[setdiff(FullList,AllParm)])
    Set2Default <- setdiff(FullList,AllParm)
  }

  # define InputPROSAIL # 1 default value
  if (length(Set2Default)>0){
    for (i in Set2Default){
      InputPROSAIL[[i]] = array(Default[i],dim = c(nbSamples,1))
    }
  }

  # define InputPROSAIL # 2 Set parameters
  if (length(ParmSet)>0){
    for (i in 1:length(ParmSet)){
      Sel <- which(InVar==names(ParmSet)[i])
      InputPROSAIL[[Sel]] = array(ParmSet[i],dim = c(nbSamples,1))
    }
  }

  # define InputPROSAIL # 3 random parameters
  for (i in 1:length(minval)){
    Sel <- names(minval)[i]
    # if uniform distribution
    if (TypeDistrib[[Sel]] == 'Uniform'){
      InputPROSAIL[[Sel]] <- array(runif(nbSamples,min = minval[1,i],max=maxval[1,i]),dim = c(nbSamples,1))
    # if Gaussian distribution
    } else if (TypeDistrib[[Sel]] == 'Gaussian'){
      InputPROSAIL[[Sel]] <- array(rnorm(2*nbSamples,mean = Mean[[Sel]],sd = Std[[Sel]]))
      Elim <- which(InputPROSAIL[[Sel]]<minval[[Sel]] | InputPROSAIL[[Sel]]>maxval[[Sel]])
      if (length(Elim)>0){
        InputPROSAIL[[Sel]] <- InputPROSAIL[[Sel]][-Elim]
      }
      if (length(InputPROSAIL[[Sel]])>=nbSamples){
        InputPROSAIL[[Sel]] <- InputPROSAIL[[Sel]][1:nbSamples]
      } else {
        # repeat to get correct number of samples
        repnb <- ceiling(nbSamples/length(InputPROSAIL[[Sel]]))
        InputPROSAIL[[Sel]] <- rep(InputPROSAIL[[Sel]],repnb)[1:nbSamples]
      }
      InputPROSAIL[[Sel]] <- array(InputPROSAIL[[Sel]],dim = c(nbSamples,1))
    }
  }
  if (Force4LowLAI==T){
    MaxLAI <- 2
    LowLAI <- which(InputPROSAIL$lai<MaxLAI)
    if (length(LowLAI)>0){
      for (BP in AllParm){
        if (!is.na(match(BP,c('CHL','CAR','EWT','LMA','ANT','PROT','CBC','N','BROWN')))){
          InputPROSAIL[[BP]][LowLAI] <- minval[[BP]]  + 0.25*(InputPROSAIL[[BP]][LowLAI]-minval[[BP]])
                                                      + InputPROSAIL$lai[LowLAI]*(0.75/MaxLAI)*(InputPROSAIL[[BP]][LowLAI]-minval[[BP]])
        }
      }
    }
  }
  return(InputPROSAIL)
}

#' This function generates distribution of biophysical parameters used as input parameters in PRO4SAIL
#'
#' @param minval list. Defines the minimum value to be set for a list of parameters randomly produced
#' @param maxval list. Defines the maximum value to be set for a list of parameters randomly produced
#' @param ParmSet list.Defines the parameters to be set to a given value
#' @param nbSamples numeric. Number of samples to be generated
#' @param TypeDistrib list. specify if uniform or Gaussian distribution to be applied. default = Uniform
#' @param Mean list. mean value for parameters with Gaussian distribution
#' @param Std list. standard deviation for parameters with Gaussian distribution
#' @param Force4LowLAI boolean. set to TRUE if needed to force distribution of leaf chemistry to low values for low LAI
#'
#' @return InputPROSAIL list. list of nbSamples input parameters for PRO4SAIL simulation
#' @importFrom stats runif rnorm sd
#' @export
get_distribution_input_prosail2 <- function(minval,maxval,ParmSet,nbSamples,
                                            TypeDistrib = data.frame('CHL'='Uniform',
                                                                     'CAR'='Uniform',
                                                                     'EWT' = 'Uniform',
                                                                     'ANT' = 'Uniform',
                                                                     'LMA' = 'Uniform',
                                                                     'PROT' = 'Uniform',
                                                                     'CBC' = 'Uniform',
                                                                     'N' = 'Uniform',
                                                                     'BROWN'='Uniform',
                                                                     'psoil' = 'Uniform',
                                                                     'LIDFa' = 'Uniform',
                                                                     'LIDFb' = 'Uniform',
                                                                     'lai' = 'Uniform',
                                                                     'q'='Uniform',
                                                                     'tto' = 'Uniform',
                                                                     'tts' = 'Uniform',
                                                                     'psi' = 'Uniform',
                                                                     'alpha' = 'Uniform',
                                                                     'fraction_brown' = 'Uniform',
                                                                     'diss' = 'Uniform',
                                                                     'Cv' = 'Uniform',
                                                                     'Zeta' = 'Uniform'),
                                            Mean = NULL,Std = NULL,Force4LowLAI=TRUE){

  # define all input parameters from PROSAIL
  InVar <- c('CHL','CAR','ANT','BROWN','EWT','LMA',
             'PROT','CBC','N','alpha','LIDFa','LIDFb',
             'lai','q','tts','tto','psi','psoil','TypeLidf',
             'fraction_brown','diss','Cv','Zeta')

  InputPROSAIL <- list('CHL'=c(),'CAR'=c(),'ANT'=c(),'BROWN'=c(),'EWT'=c(),
                       'LMA'=c(),'PROT'=c(),'CBC'=c(),'N'=c(),'alpha'=c(),
                       'LIDFa'=c(),'LIDFb'=c(),'lai'=c(),'q'=c(),
                       'tts'=c(),'tto'=c(),'psi'=c(),'psoil'=c(),'TypeLidf'=c(),
                       'fraction_brown'=c(),'diss'=c(),'Cv'=c(),'Zeta'=c())
  Default <- data.frame('CHL'=0,'CAR'=0,'ANT'=0,'BROWN'=0,'EWT'=0,
                        'LMA'=0,'PROT'=0,'CBC'=0,'N'=1.5,'alpha'=40,
                        'LIDFa'=0,'LIDFb'=0,'lai'=2,'q'=0,
                        'tts'=0,'tto'=0,'psi'=0,'psoil'=1,'TypeLidf'=2,
                        'fraction_brown'=0.5,'diss'=0.5,'Cv'=1,'Zeta'=1)

  # which input parameters should be randomly sampled?
  ParmRand <- which(is.element(InVar,names(minval))==TRUE)
  # which input parameters should be set to a constant value?
  ParmCte <- which(is.element(InVar,names(ParmSet))==TRUE)
  # if some parameters are defined to be sampled randomly and set to constant value
  if (length(intersect(ParmRand,ParmCte))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(InVar[intersect(ParmRand,ParmCte)])
    message('will be randomly sampled based on the min and max values defined by minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  AllParm <- c(ParmRand,ParmCte)
  FullList <- seq(1,length(InVar),by=1)
  Set2Default <- c()
  if (length(setdiff(FullList,AllParm))>0){
    message('WARNING: some elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(InVar[setdiff(FullList,AllParm)])
    message('These parameters will be set to their default value')
    print(Default[setdiff(FullList,AllParm)])
    Set2Default <- setdiff(FullList,AllParm)
  }

  # define InputPROSAIL # 1 default value
  if (length(Set2Default)>0){
    for (i in Set2Default){
      InputPROSAIL[[i]] = array(Default[i],dim = c(nbSamples,1))
    }
  }

  # define InputPROSAIL # 2 Set parameters
  if (length(ParmSet)>0){
    for (i in 1:length(ParmSet)){
      Sel <- which(InVar==names(ParmSet)[i])
      InputPROSAIL[[Sel]] = array(ParmSet[i],dim = c(nbSamples,1))
    }
  }

  # define InputPROSAIL # 3 random parameters
  for (i in 1:length(minval)){
    Sel <- names(minval)[i]
    # if uniform distribution
    if (TypeDistrib[[Sel]] == 'Uniform'){
      InputPROSAIL[[Sel]] <- array(runif(nbSamples,min = minval[1,i],max=maxval[1,i]),dim = c(nbSamples,1))
      # if Gaussian distribution
    } else if (TypeDistrib[[Sel]] == 'Gaussian'){
      InputPROSAIL[[Sel]] <- array(rnorm(2*nbSamples,mean = Mean[[Sel]],sd = Std[[Sel]]))
      Elim <- which(InputPROSAIL[[Sel]]<minval[[Sel]] | InputPROSAIL[[Sel]]>maxval[[Sel]])
      if (length(Elim)>0){
        InputPROSAIL[[Sel]] <- InputPROSAIL[[Sel]][-Elim]
      }
      if (length(InputPROSAIL[[Sel]])>=nbSamples){
        InputPROSAIL[[Sel]] <- InputPROSAIL[[Sel]][1:nbSamples]
      } else {
        # repeat to get correct number of samples
        repnb <- ceiling(nbSamples/length(InputPROSAIL[[Sel]]))
        InputPROSAIL[[Sel]] <- rep(InputPROSAIL[[Sel]],repnb)[1:nbSamples]
      }
      InputPROSAIL[[Sel]] <- array(InputPROSAIL[[Sel]],dim = c(nbSamples,1))
    }
  }

  if (Force4LowLAI==T){
    MaxLAI <- 2
    LowLAI <- which(InputPROSAIL$lai<MaxLAI)
    if (length(LowLAI)>0){
      for (BP in InVar){
        if (!is.na(match(BP,c('CHL','CAR','EWT','LMA','ANT','PROT','CBC','N','BROWN')))){
          if (length(unique(InputPROSAIL[[BP]]))>1){
            InputPROSAIL[[BP]][LowLAI] <- minval[[BP]]  + 0.20*(InputPROSAIL[[BP]][LowLAI]-minval[[BP]]) +
                                                          InputPROSAIL$lai[LowLAI]*(0.80/MaxLAI)*(InputPROSAIL[[BP]][LowLAI]-minval[[BP]])
          }
        }
      }
    }
  }
  return(InputPROSAIL)
}

#' This function generates a LUT of BRF based on a table of input variables for PRO4SAIL model
#'
#' @param InputPROSAIL list. PROSAIL input variables
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for clear conditions
#' @param BandNames character. Name of the spectral bands of the sensor
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param BrownVegetation list. Defines optical properties for brown vegetation, if not NULL
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#'
#' @return BRF_LUT numeric. matrix of BRF corresponding to InputPROSAIL
#' @importFrom progress progress_bar
#' @export

Generate_LUT_BRF <- function(InputPROSAIL, SpecPROSPECT, SpecSOIL, SpecATM, BandNames = NULL,
                             SAILversion='4SAIL', BrownVegetation = NULL){

  nbSamples <- length(InputPROSAIL[[1]])
  BRF <- list()
  Split <- round(nbSamples/10)
  pb <- progress_bar$new(
    format = "Generate LUT [:bar] :percent in :elapsed",
    total = 10, clear = FALSE, width= 100)
  for (i in 1:nbSamples){
    if (i%%Split==0 & nbSamples>100){
        pb$tick()
    }
    rsoil <- InputPROSAIL$psoil[[i]]*SpecSOIL$Dry_Soil+(1-InputPROSAIL$psoil[[i]])*SpecSOIL$Wet_Soil
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,CHL = InputPROSAIL$CHL[[i]], CAR = InputPROSAIL$CAR[[i]],
                          ANT = InputPROSAIL$ANT[[i]], EWT = InputPROSAIL$EWT[[i]], LMA = InputPROSAIL$LMA[[i]],
                          PROT = InputPROSAIL$PROT[[i]], CBC = InputPROSAIL$CBC[[i]], BROWN = InputPROSAIL$BROWN[[i]],
                          N = InputPROSAIL$N[[i]],
                          TypeLidf = InputPROSAIL$TypeLidf[[i]],LIDFa = InputPROSAIL$LIDFa[[i]],LIDFb = InputPROSAIL$LIDFb[[i]],
                          lai = InputPROSAIL$lai[[i]],q = InputPROSAIL$q[[i]],
                          tts = InputPROSAIL$tts[[i]],tto = InputPROSAIL$tto[[i]],psi = InputPROSAIL$psi[[i]],rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,CHL = InputPROSAIL$CHL[[i]], CAR = InputPROSAIL$CAR[[i]],
                          ANT = InputPROSAIL$ANT[[i]], EWT = InputPROSAIL$EWT[[i]], LMA = InputPROSAIL$LMA[[i]],
                          PROT = InputPROSAIL$PROT[[i]], CBC = InputPROSAIL$CBC[[i]], BROWN = InputPROSAIL$BROWN[[i]],
                          N = InputPROSAIL$N[[i]],
                          TypeLidf = InputPROSAIL$TypeLidf[[i]],LIDFa = InputPROSAIL$LIDFa[[i]],LIDFb = InputPROSAIL$LIDFb[[i]],
                          lai = InputPROSAIL$lai[[i]],q = InputPROSAIL$q[[i]],
                          tts = InputPROSAIL$tts[[i]],tto = InputPROSAIL$tto[[i]],psi = InputPROSAIL$psi[[i]],rsoil = rsoil,
                          SAILversion='4SAIL2',
                          fraction_brown = InputPROSAIL$fraction_brown[[i]], diss = InputPROSAIL$diss[[i]], Cv = InputPROSAIL$Cv[[i]],Zeta = InputPROSAIL$Zeta[[i]],
                          BrownVegetation = BrownVegetation)
    }
    # Computes bidirectional reflectance factor based on outputs from PROSAIL and sun position
    BRF[[i]] <-Compute_BRF(RefSAIL$rdot,RefSAIL$rsot,InputPROSAIL$tts[[i]],SpecATM)
  }
  BRF <- do.call(cbind,BRF)
  row.names(BRF) <- BandNames
  return(BRF)
}

#' This function applied noise on a matrix
#'
#' @param LUT numeric. Matrix including data to add noise to
#' @param NoiseLevel numeric. value of the normal noise proportional to LUT to apply on LUT
#' @param NoiseType character.
#' - relative: noise proportional to actual value to add noise to
#' - absolute: noise not proportional to actual value to add noise to
#'
#' @return LUT_Noise numeric. Matrix including data with added noise
#' @export

Apply_Noise_LUT <- function(LUT,NoiseLevel,NoiseType = 'relative'){
  nbFeatures <- nrow(LUT)
  nbSamples <- ncol(LUT)
  if (NoiseType == 'relative'){
    LUT_Noise <- LUT + LUT*matrix(rnorm(nbFeatures*nbSamples,0,NoiseLevel),nrow = nbFeatures)
  } else if (NoiseType == 'absolute'){
    LUT_Noise <- LUT + matrix(rnorm(nbFeatures*nbSamples,0,NoiseLevel),nrow = nbFeatures)
  }
  return(LUT_Noise)
}



