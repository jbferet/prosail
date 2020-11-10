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

#' This function generates distribution of biophysical parameters used as input parameters in PRO4SAIL
#'
#' @param minval list. Defines the minimum value to be set for a list of parameters randomly produced
#' @param maxval list. Defines the maximum value to be set for a list of parameters randomly produced
#' @param ParmSet list.Defines the parameters to be set to a given value
#' @param nbSamples numeric. Number of samples to be generated
#' @param TypeDistrib list. specify if uniform or Gaussian distribution to be applied. default = Uniform
#' @param Mean list. mean value for parameters with Gaussian distribution
#' @param Std list. standard deviation for parameters with Gaussian distribution
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
                                           Mean = NULL,Std = NULL){

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
  return(InputPROSAIL)
}


#' This function generates a LUT of BRF based on a table of input variables for PRO4SAIL model
#'
#' @param InputPROSAIL list. PROSAIL input variables
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for clear conditions
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param BrownVegetation list. Defines optical properties for brown vegetation, if not NULL
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#'
#' @return BRF_LUT numeric. matrix of BRF corresponding to InputPROSAIL
#' @importFrom progress progress_bar
#' @export

Generate_LUT_BRF <- function(InputPROSAIL,SpecPROSPECT,SpecSOIL,SpecATM,SAILversion='4SAIL',BrownVegetation = NULL){

  nbSamples <- length(InputPROSAIL[[1]])
  BRF <- list()
  Split <- round(nbSamples/10)
  pb <- progress_bar$new(
    format = "Generate LUT [:bar] :percent in :elapsed",
    total = 10, clear = FALSE, width= 100)
  for (i in 1:nbSamples){
    if (i%%Split==0){
      pb$tick()
      Sys.sleep(1 / 10)
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
  return(BRF)
}

#' This function applied noise on a matrix
#'
#' @param LUT numeric. Matrix including data to add noise to
#' @param NoiseLevel numeric. value of teh normal noise proportional to LUT to apply on LUT
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
  } else if (NoiseType == 'relative'){
    LUT_Noise <- LUT + matrix(rnorm(nbFeatures*nbSamples,0,NoiseLevel),nrow = nbFeatures)
  }
  return(LUT_Noise)
}



