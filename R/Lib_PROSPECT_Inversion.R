# ==============================================================================
# prospect
# Lib_PROSPECT_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to PROSPECT inversion
# ==============================================================================

#' Performs PROSPECT inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param Parms2Estimate  list. Parameters to estimate (can be 'ALL')
#' @param ParmSet  list. Parameters set to a fixed value by user
#' @param PROSPECT_version  character. '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom pracma fmincon
#' @export
Invert_PROSPECT  <- function(SpecPROSPECT,Refl = NULL,Tran = NULL,
                             Parms2Estimate = 'ALL',ParmSet = NULL,
                             PROSPECT_version = 'D',MeritFunction = 'Merit_RMSE'){

  # define PROSPECT input parameters
  InPROSPECT  = data.frame('CHL'=0,'CAR'=0,'ANT'=0,'BROWN'=0,'EWT'=0,
                           'LMA'=0,'PROT'=0,'CBC'=0,'N'=0,'alpha'=40)
  # update InPROSPECT according to ParmSet
  if (!is.null(ParmSet)){
    for (i in 1:length(names(ParmSet))){
      InPROSPECT[which(names(InPROSPECT)%in%names(ParmSet)[i])] = ParmSet[i]
    }
  }
  # define init value and lower/upper boundaries for inversion
  #               CHL   CAR   ANT   BROWN   EWT   LMA   PROT  CBC   N     alpha
  xinit_All   = c(40,   10,   0.1,  0.01,   0.01, 0.008,0.001,0.008,    1.5,  40)
  lb_All      = c(1e-4, 1e-4, 0,    0,      1e-7, 1e-7, 0.0,  0.00,     0.5,  10)
  ub_All      = c(150,  25,   20.0, 1.0,    0.08, 0.04, 0.005,0.04,     3.0,  90)
  # # update init value and lower/upper boundaries for inversion based on user values
  # if (!is.null(xinit_user)){
  # }

  # set parameters to user value defined in ParmSet
  Vars2Estimate = c()
  if ('ALL'%in%Parms2Estimate){
    if (PROSPECT_version == '5'){
      Vars2Estimate = c(1,2,5,6,9)
    } else if (PROSPECT_version == '5B'){
      Vars2Estimate = c(1,2,4,5,6,9)
    } else if (PROSPECT_version == 'D'){
      Vars2Estimate = c(1,2,3,5,6,9)
    } else if (PROSPECT_version == 'DB'){
      Vars2Estimate = c(1,2,3,4,5,6,9)
    } else if (PROSPECT_version == 'PRO'){
      Vars2Estimate = c(1,2,3,5,7,8,9)
    } else if (PROSPECT_version == 'PROB'){
      Vars2Estimate = c(1,2,3,4,5,7,8,9)
    }
  } else {
    if ('CHL'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,1)}
    if ('CAR'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,2)}
    if ('ANT'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,3)}
    if ('BROWN'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,4)}
    if ('EWT'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,5)}
    if ('LMA'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,6)}
    if ('PROT'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,7)}
    if ('CBC'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,8)}
    if ('N'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,9)}
    if (!'N'%in%Parms2Estimate & !'N'%in%names(ParmSet)){Vars2Estimate = c(Vars2Estimate,9)}
    if ('alpha'%in%Parms2Estimate){Vars2Estimate = c(Vars2Estimate,10)}
  }

  # update init value and lower/upper boundaries for inversion based on Vars2Estimate
  xinit = xinit_All[Vars2Estimate]
  lb    = lb_All[Vars2Estimate]
  ub    = ub_All[Vars2Estimate]
  res   = fmincon(x0 = xinit, fn = MeritFunction, gr = NULL,
                  SpecPROSPECT=SpecPROSPECT,Refl=Refl,Tran=Tran,
                  Input_PROSPECT = InPROSPECT,WhichVars2Estimate=Vars2Estimate,
                  method = "SQP",A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                  lb = lb, ub = ub, hin = NULL, heq = NULL,tol = 1e-07,
                  maxfeval = 2000, maxiter = 1000)

  OutPROSPECT  = InPROSPECT
  OutPROSPECT[Vars2Estimate]= res$par
  return(OutPROSPECT)
}

#' Merit function for PROSPECT-D and estimation of EWT or LMA
#'
#' @param xinit numeric. Vector of input variables to estimate
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param Input_PROSPECT dataframe. set of PROSPECT input variables
#' @param WhichVars2Estimate  numeric. location of variables from Input_PROSPECT
#' to be estimated through inversion
#'
#' @return fc estimates of the parameters
#' @export
Merit_RMSE <- function(xinit,SpecPROSPECT,Refl,Tran,Input_PROSPECT,WhichVars2Estimate) {

  xinit[xinit<0] = 0
  Input_PROSPECT[WhichVars2Estimate] = xinit
  RT = PROSPECT(SpecPROSPECT = SpecPROSPECT,Input_PROSPECT = Input_PROSPECT)
  fc = CostVal_RMSE(RT,Refl,Tran)
  return(fc)
}


#' Value of the cost criterion to minimize during PROSPECT inversion
#' @param RT  list. Simulated reflectance and transmittance
#' @param Refl  numeric. Reflectance on which PROSPECT ins inverted
#' @param Tran  numeric. Transmittance on which PROSPECT ins inverted
#'
#' @return fc sum of squared difference between simulated and measured leaf optical properties
#' @export
CostVal_RMSE  <- function(RT,Refl,Tran) {

  if (is.null(Tran)){
    fc = sqrt(sum((Refl-RT$Reflectance)**2)/length(RT$Reflectance))
  } else if (is.null(Refl)){
    fc = sqrt(sum((Tran-RT$Transmittance)**2)/length(RT$Transmittance))
  } else {
    fc=sqrt(sum((Refl-RT$Reflectance)**2)/length(RT$Reflectance)+sum((Tran-RT$Transmittance)**2)/length(RT$Transmittance))
  }
  return(fc)
}

#' This function adapts SpecPROSPECT accordingly to experimental data
#' or to a spectral domain defined by UserDomain
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param lambda  numeric. Spectral bands corresponding to experimental data
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param UserDomain  numeric. Lower and upper bounds for domain of interest (optional)
#'
#' @return res list including spectral properties at the new resolution
#' @importFrom utils tail
#' @export
FitSpectralData <- function(SpecPROSPECT,lambda,Refl=NULL,Tran=NULL,UserDomain=NULL) {

  LowerPROSPECT = SpecPROSPECT$lambda[1]
  UpperPROSPECT = tail(SpecPROSPECT$lambda, n=1)
  LowerLOP      = lambda[1]
  UpperLOP      = tail(lambda, n=1)
  # if only need to fit PROSPECT data with spctral dat aprovided by user
  if (is.null(UserDomain)){
    if (LowerPROSPECT>LowerLOP){
      warning('________________________ WARNING _______________________')
      warning('User spectal data will be shrinked to start at the same ')
      warning('         Spectral band as SpecPROSPECT, which is        ')
      print(LowerPROSPECT)
      LowerBand_LOP = which(abs(lambda-LowerPROSPECT)==min(abs(lambda-LowerPROSPECT)))
      LowerBand_Spec= 1
    } else  if (LowerPROSPECT<LowerLOP){
      LowerBand_Spec  = which(abs(SpecPROSPECT$lambda-LowerLOP)==min(abs(SpecPROSPECT$lambda-LowerLOP)))
      LowerBand_LOP = 1
    } else  if (LowerPROSPECT==LowerLOP){
      LowerBand_Spec= 1
      LowerBand_LOP = 1
    }
    if (UpperPROSPECT<UpperLOP){
      warning('________________________ WARNING _______________________')
      warning('  User spectal data will be shrinked to end at the same ')
      warning('         Spectral band as SpecPROSPECT, which is        ')
      print(UpperPROSPECT)
      UpperBand_LOP = which(abs(lambda-UpperPROSPECT)==min(abs(lambda-UpperPROSPECT)))
      UpperBand_Spec= length(SpecPROSPECT$lambda)
    } else  if (UpperPROSPECT>UpperLOP){
      UpperBand_Spec  = which(abs(SpecPROSPECT$lambda-UpperLOP)==min(abs(SpecPROSPECT$lambda-UpperLOP)))
      UpperBand_LOP = length(lambda)
    } else  if (UpperPROSPECT==UpperLOP){
      UpperBand_Spec= length(SpecPROSPECT$lambda)
      UpperBand_LOP = length(lambda)
    }
    # if user specifies a spectral domain which is different from PROSPECT and user data
  } else if (!is.null(UserDomain)){
    LowerUser     = UserDomain[1]
    UpperUser     = UserDomain[2]
    if (LowerLOP>LowerUser | UpperLOP<UpperUser | LowerPROSPECT>LowerUser | UpperPROSPECT<UpperUser){
      if (LowerPROSPECT>LowerUser | UpperPROSPECT<UpperUser){
        warning('________________________ WARNING _______________________')
        warning('  The spectral domain defined in UserDomain provided as ')
        warning(' input in function FitSpectralData does not match with  ')
        warning('       the spectral domain covered by PROSPECT          ')
        warning('                                                        ')
        warning('                 PLEASE ADJUST UserDomain               ')
        warning('                                                        ')
        stop()
      }
      if (LowerLOP>LowerUser | UpperLOP<UpperUser){
        warning('________________________ WARNING _______________________')
        warning('  The spectral domain defined in UserDomain provided as ')
        warning(' input in function FitSpectralData does not match with  ')
        warning('       the spectral domain covered by user data         ')
        warning('                                                        ')
        warning('                 PLEASE ADJUST UserDomain               ')
        warning('                                                        ')
        stop()
      }
    } else {
      if (LowerLOP<=LowerUser){
        LowerBand_LOP  = which(abs(lambda-LowerUser)==min(abs(lambda-LowerUser)))
      }
      if (UpperLOP>=UpperUser){
        UpperBand_LOP  = which(abs(lambda-UpperUser)==min(abs(lambda-UpperUser)))
      }
      if (LowerPROSPECT<=LowerUser){
        LowerBand_Spec  = which(abs(SpecPROSPECT$lambda-LowerUser)==min(abs(SpecPROSPECT$lambda-LowerUser)))
      }
      if (UpperPROSPECT>=UpperUser){
        UpperBand_Spec  = which(abs(SpecPROSPECT$lambda-UpperUser)==min(abs(SpecPROSPECT$lambda-UpperUser)))
      }
    }
  }
  SubSpecPROSPECT = SpecPROSPECT[LowerBand_Spec:UpperBand_Spec,]
  Sublambda       = lambda[LowerBand_LOP:UpperBand_LOP]
  if (!length(Sublambda)==length(SubSpecPROSPECT$lambda)){
    warning('______________________ WARNING _____________________')
    warning('       PROSPECT expects 1nm spectal sampling        ')
    warning('The data provided as input shows unexpected sampling')
    warning('    Please prepare your data accordingly before     ')
    warning('             running PROSPECT inversion             ')
    warning('                The process will stop               ')
    stop()
  }
  SubRefl = SubTran = NULL
  if (!is.null(Refl)){
    SubRefl = Refl[LowerBand_LOP:UpperBand_LOP,]
  }
  if (!is.null(Tran)){
    SubTran = Tran[LowerBand_LOP:UpperBand_LOP,]
  }
  res       = list('SpecPROSPECT'=SubSpecPROSPECT,'lambda'=Sublambda,'Refl'=SubRefl,'Tran'=SubTran)
  return(res)
}


#' This function defines a regression model to estimate N from R only or T only
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#'
#' @return Nprior vector corresponding to teh prior estimation of N based on R only or T only
#' @importFrom stats lm runif
#' @export
Get_Nprior <- function(SpecPROSPECT,lambda,Refl=NULL,Tran=NULL) {

  # definition of the optimal spectral band based on data available
  OptWL_R = OptWL_T = list()
  OptWL_R$NIR     = 800
  OptWL_R$SWIR    = 1131
  OptWL_T$NIR     = 753
  OptWL_T$SWIR    = 1121
  # if prior information based on Reflectance
  if (is.null(Tran)){
    if (OptWL_R$SWIR%in%SpecPROSPECT$lambda){
      OptWL = OptWL_R$SWIR
    } else if (!OptWL_R$SWIR%in%SpecPROSPECT$lambda & OptWL_R$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The optimal prior estimation of N using Reflectance only')
      warning('requires information at 1131nm.')
      warning('The reflectance does not include this spectral band')
      warning('Using reflectance at 800 nm instead')
      OptWL = OptWL_R$NIR
    } else if (!OptWL_R$SWIR%in%SpecPROSPECT$lambda & !OptWL_R$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The spectral information of the reflectance provided here')
      warning('does not contain the spectral bands required to estimate')
      warning('prior information about N.')
      warning('The proces will stop')
      stop()
    }
  } else if (is.null(Refl)){
    if (OptWL_T$SWIR%in%SpecPROSPECT$lambda){
      OptWL = OptWL_T$SWIR
    } else if (!OptWL_T$SWIR%in%SpecPROSPECT$lambda & OptWL_T$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The optimal prior estimation of N using Transmittance only')
      warning('requires information at 1121nm.')
      warning('The Transmittance does not include this spectral band')
      warning('Using Transmittance at 753 nm instead')
      OptWL = OptWL_T$NIR
    } else if (!OptWL_T$SWIR%in%SpecPROSPECT$lambda & !OptWL_T$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The spectral information of the Transmittance provided here')
      warning('does not contain the spectral bands required to estimate')
      warning('prior information about N.')
      warning('The proces will stop')
      stop()
    }
  }

  # get the subdomain corresponding to OptWL
  SubData = FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = c(OptWL,OptWL))
  SubSpecPROSPECT = SubData$SpecPROSPECT
  Sublambda       = SubData$lambda
  SubRefl         = SubData$Refl
  subTran         = SubData$Tran

  # create a LUT using only the spectral band of interest
  nbSim   = 1000
  CHL     = 0.5+100*runif(nbSim)
  CAR     = 0.5+20*runif(nbSim)
  EWT     = 0.001+0.02*runif(nbSim)
  LMA     = 0.001+0.01*runif(nbSim)
  N       = 1+1.5*runif(nbSim)
  Input_PROSPECT = data.frame('CHL'=CHL,'CAR'=CAR,'EWT'=EWT,'LMA'=LMA,'N'=N)
  LUT   = PROSPECT_LUT(SubSpecPROSPECT,Input_PROSPECT)
  # fit a linear model between
  if (is.null(Tran)){
    Ratio       = LUT$Reflectance/(1-LUT$Reflectance)
    Ratio_Meas  = SubRefl/(1-SubRefl)
  } else if (is.null(Refl)){
    Ratio       = (1-LUT$Transmittance)/LUT$Transmittance
    Ratio_Meas  = (1-subTran)/subTran
  }
  N_Model   = lm(matrix(LUT$Input_PROSPECT$N) ~ matrix(Ratio))
  NpriorMOD = N_Model$coefficients[2]*Ratio+N_Model$coefficients[1]
  Nprior    = N_Model$coefficients[2]*Ratio_Meas+N_Model$coefficients[1]
  return(Nprior)
}
