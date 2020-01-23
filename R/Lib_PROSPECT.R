# ==============================================================================
# prospect
# Lib_PROSPECT.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to PROSPECT simulation
# ==============================================================================

#' core function running PROSPECT
#  this code includes numerical optimizations proosed in the FLUSPECT code
#  Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans,
#  Date: 2007
#  Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
#'
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Input_PROSPECT dataframe. full list of PROSPECT input parameters. Default = Set to NULL
#' @param N numeric. Leaf structure parameter
#' @param CHL numeric. Chlorophyll content (microg.cm-2)
#' @param CAR numeric. Carotenoid content (microg.cm-2)
#' @param ANT numeric. Anthocyain content (microg.cm-2)
#' @param BROWN numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param PROT numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf (simulation of roughness)
#'
#' @return leaf directional-hemisphrical reflectance and transmittance
#' @importFrom expint expint
#' @export
PROSPECT  <- function(SpecPROSPECT,Input_PROSPECT=NULL,N = 1.5,CHL = 40.0,
                      CAR = 8.0,ANT = 0.0,BROWN = 0.0,EWT = 0.01,
                      LMA = 0.008,PROT = 0.0,CBC = 0.0,alpha = 40.0){

  # if PROSPECT values are provided as individual parametrs
  if (is.null(Input_PROSPECT)){
    Input_PROSPECT = data.frame('CHL'= CHL, 'CAR'= CAR, 'ANT'=ANT, 'BROWN'= BROWN, 'EWT'=EWT,
                                'LMA'=LMA, 'PROT'= PROT, 'CBC'= CBC, 'N'=N, 'alpha'= alpha)
  }
  if (is.null(Input_PROSPECT$alpha)){
    Input_PROSPECT$alpha = alpha
  }

  if (!is.null(Input_PROSPECT$PROT) | !is.null(Input_PROSPECT$CBC)){
    if (Input_PROSPECT$LMA>0 & (Input_PROSPECT$PROT >0 | Input_PROSPECT$CBC >0)){
      message('PROT and/or CBC are not set to 0')
      message('LMA is not set to 0 neither, which is physically incorrect')
      message('(LMA = PROT + CBC)')
      message('We assume that PROSPECT-PRO was called and set LMA to 0')
      message('Please correct input parameters LMA, PROT and/or CBC if needed')
      Input_PROSPECT$LMA <- 0
    }
  }

  Kall <-   (Input_PROSPECT$CHL*SpecPROSPECT$SAC_CHL + Input_PROSPECT$CAR*SpecPROSPECT$SAC_CAR+
               Input_PROSPECT$ANT*SpecPROSPECT$SAC_ANT + Input_PROSPECT$BROWN*SpecPROSPECT$SAC_BROWN+
               Input_PROSPECT$EWT*SpecPROSPECT$SAC_EWT + Input_PROSPECT$LMA*SpecPROSPECT$SAC_LMA+
               Input_PROSPECT$PROT*SpecPROSPECT$SAC_PROT+ Input_PROSPECT$CBC*SpecPROSPECT$SAC_CBC)/Input_PROSPECT$N

  j           = which(Kall>0)                # Non-conservative scattering (normal case)
  t1          = (1-Kall)*exp(-Kall)
  t2          = (Kall*Kall)*expint(Kall)
  tau         = matrix(1,ncol = 1,nrow = length(t1))
  tau[j]      = t1[j]+t2[j]

  # ***********************************************************************
  # reflectance and transmittance of one layer
  # ***********************************************************************
  # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
  # Interaction of isotropic ligth with a compact plant leaf, J. Opt.
  # Soc. Am., 59(10):1376-1379.
  # ***********************************************************************
  # reflectivity and transmissivity at the interface
  # ***********************************************************************
  talf        = calctav(Input_PROSPECT$alpha,SpecPROSPECT$nrefrac)
  ralf        = 1-talf
  t12         = calctav(90,SpecPROSPECT$nrefrac)
  r12         = 1-t12
  t21         = t12/(SpecPROSPECT$nrefrac**2)
  r21         = 1-t21

  # top surface side
  denom       = 1-(r21*r21*(tau**2))
  Ta          = (talf*tau*t21)/denom
  Ra          = ralf+(r21*tau*Ta)
  # bottom surface side
  t           = t12*tau*t21/denom
  r           = r12+(r21*tau*t)

  # ***********************************************************************
  # reflectance and transmittance of N layers
  # Stokes equations to compute properties of next N-1 layers (N real)
  # Normal case
  # ***********************************************************************
  # Stokes G.G. (1862), On the intensity of the light reflected from
  # or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
  # 11:545-556.
  # ***********************************************************************
  D           = sqrt((1+r+t)*(1+r-t)*(1-r+t)*(1-r-t))
  rq          = r**2
  tq          = t**2
  a           = (1+rq-tq+D)/(2*r)
  b           = (1-rq+tq+D)/(2*t)

  bNm1        = b**(Input_PROSPECT$N-1)
  bN2         = bNm1**2
  a2          = a**2
  denom       = a2*bN2-1
  Rsub        = a*(bN2-1)/denom
  Tsub        = bNm1*(a2-1)/denom

  #	Case of zero absorption
  j           = which(r+t >= 1)
  Tsub[j]     = t[j]/(t[j]+(1-t[j])*(N-1))
  Rsub[j]	    = 1-Tsub[j]

  # Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
  denom       = 1-Rsub*r
  tran        = Ta*Tsub/denom
  refl        = Ra+(Ta*Rsub*t)/denom
  my_list <- list("Reflectance" = refl,"Transmittance" =tran,"wvl"=SpecPROSPECT$lambda)
  return(my_list)
}

#' computation of transmissivity of a dielectric plane surface,
#' averaged over all directions of incidence and over all polarizations.
#'
#' @param alpha numeric. Maximum incidence angle defining the solid angle of incident light
#' @param nr numeric. refractive index
#'
#' @return numeric. Transmissivity of a dielectric plane surface
#' @export
calctav  <- function(alpha,nr) {
  # Stern F. (1964), Transmission of isotropic radiation across an
  # interface between two dielectrics, Appl. Opt., 3(1):111-113.
  # Allen W.A. (1973), Transmission of isotropic light across a
  # dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
  # 63(6):664-666.
  # ***********************************************************************

  rd  = pi/180
  n2  = nr**2
  np  = n2+1
  nm  = n2-1
  a   = (nr+1)*(nr+1)/2
  k   = -(n2-1)*(n2-1)/4
  sa  = sin(alpha*rd)

  b2  = (sa**2)-(np/2)
  if (alpha==90){
    b1  = 0*b2
  } else {
    b1  = sqrt((b2**2)+k)
  }
  b   = b1-b2
  b3  = b**3
  a3  = a**3
  ts  = ((k**2)/(6*b3)+(k/b)-b/2)-((k**2)/(6*a3)+(k/a)-(a/2))

  tp1 = -2*n2*(b-a)/(np**2)
  tp2 = -2*n2*np*log(b/a)/(nm**2)
  tp3 = n2*((1/b)-(1/a))/2
  tp4 = 16*n2**2*((n2**2)+1)*log(((2*np*b)-(nm**2))/((2*np*a)-(nm**2)))/((np**3)*(nm**2))
  tp5 = 16*(n2**3)*(1/((2*np*b)-(nm**2))-(1/(2*np*a-(nm**2))))/(np**3)
  tp  = tp1+tp2+tp3+tp4+tp5
  tav = (ts+tp)/(2*(sa**2))
  return(tav)
}

#' computation of a LUT of leaf optical properties using a set of
#' leaf chemical & structural parameters
#'
#' @param SpecPROSPECT dataframe. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Input_PROSPECT dataframe. full list of PROSPECT input parameters. Default = Set to NULL
#' @param N numeric. Leaf structure parameter
#' @param CHL numeric. Chlorophyll content (microg.cm-2)
#' @param CAR numeric. Carotenoid content (microg.cm-2)
#' @param ANT numeric. Anthocyain content (microg.cm-2)
#' @param BROWN numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param PROT numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf (simulation of roughness)
#'
#' @return list. LUT including leaf reflectance and transmittance
#' @export
PROSPECT_LUT  <- function(SpecPROSPECT,Input_PROSPECT=NULL,N = NULL,CHL = NULL,
                          CAR = NULL,ANT = NULL,BROWN = NULL,EWT = NULL,
                          LMA = NULL,PROT = NULL,CBC = NULL,alpha = NULL){

  # expected PROSPECT input parameters
  ExpectedParms =   data.frame('CHL'=0,'CAR'=0,'ANT'=0,'BROWN'=0,'EWT'=0,
                               'LMA'=0,'PROT'=0,'CBC'=0,'N'=1.5,'alpha'=40)
  # if parameters already provided as a list
  if (!is.null(Input_PROSPECT)){
    # identify missing elements
    Parm2Add = which(names(ExpectedParms)%in%names(Input_PROSPECT)==FALSE)
    # check if all parameters are included.
    if (length(Parm2Add)>0){
      # print warning
      Warn_MissingInput(Input_PROSPECT,ExpectedParms)
    }
  } else if (is.null(Input_PROSPECT)){
    # create a list of input parameters
    Input_PROSPECT  = list('CHL'=CHL,'CAR'=CAR,'ANT'=ANT,
                           'BROWN'=BROWN,'EWT'=EWT,'LMA'=LMA,
                           'PROT'=PROT,'CBC'=CBC,'N'=N,
                           'alpha'=alpha)
    # check if one or some input parametrs are NULL
    OK_Parm         = which((lengths(Input_PROSPECT) != 0)==TRUE)
    Parm2Add       = which((lengths(Input_PROSPECT) != 0)==FALSE)
    # if missing parameter
    if (length(Parm2Add)>0){
      Input_PROSPECT = Input_PROSPECT[-Parm2Add]
      # print warning
      Warn_MissingInput(Input_PROSPECT,ExpectedParms)
    }
  }
  # re-order missing elements the end of the list using default value
  Input_PROSPECT = Complete_Input_PROSPECT(Input_PROSPECT,Parm2Add,ExpectedParms)
  # print number of samples to be simulated
  nbSamples = length(Input_PROSPECT[[1]])
  paste('A LUT will be produced, including', as.character(nbSamples), 'samples', sep=' ')

  # run PROSPECT for nbSamples
  LUT_Refl = LUT_Tran = matrix(0,nrow = length(SpecPROSPECT$lambda),ncol =nbSamples)
  for (i in 1:nbSamples){
    LUT_tmp   = PROSPECT(SpecPROSPECT,Input_PROSPECT[i,])
    LUT_Refl[,i]  = LUT_tmp$Reflectance
    LUT_Tran[,i]  = LUT_tmp$Transmittance
  }
  LUT = list('Reflectance'=LUT_Refl, 'Transmittance'=LUT_Tran, 'Input_PROSPECT'=Input_PROSPECT)
  return(LUT)
}

#' Display warning message to inform about the parametrs default values
#' used during PROSPECT simulations
#'
#' @param Input_PROSPECT list. Input parameters sent to PROSPECT by user
#' @param ExpectedParms list. Full set of parameters expected to run PROSPECT
#'
#' @return WhichMissing missing parameters
#' @export
Warn_MissingInput  <- function(Input_PROSPECT,ExpectedParms){
  message('_________________ WARNING _________________')
  message('  The list of parameters provided for LUT  ')
  message('      in the Input_PROSPECT variable       ')
  message('does not include the full set of parameters')
  message('')
  # message('        Set of parameters included:        ')
  # print(names(Input_PROSPECT))
  # message('')
  # message('        Set of parameters expected:        ')
  # print(names(ExpectedParms))
  # message('')
  message(' The following parameters will be set to  ')
  message('   their default value as given here      ')
  WhichMissing = (ExpectedParms[which(names(ExpectedParms)%in%names(Input_PROSPECT)==FALSE)])
  print(WhichMissing)
}

#' Complete the list of PROSPECT parameters with default values
#'
#' @param Input_PROSPECT input parameters sent to PROSPECT by user
#' @param Parm2Add Parameters to be added to input parameters
#' @param ExpectedParms full set of parameters expected to run PROSPECT
#'
#' @return Input_PROSPECT
#' @export
Complete_Input_PROSPECT  <- function(Input_PROSPECT,Parm2Add,ExpectedParms){
  ii=0
  nbSamples = length(Input_PROSPECT[[1]])
  nbInputs  = length(Input_PROSPECT)
  for (i in Parm2Add){
    ii=ii+1
    nbInputs = nbInputs+1
    Input_PROSPECT[[nbInputs]] = matrix(ExpectedParms[[i]],ncol = 1,nrow = nbSamples)
    names(Input_PROSPECT)[[nbInputs]] = names(ExpectedParms)[[i]]
  }
  Input_PROSPECT  = data.frame('CHL'=matrix(Input_PROSPECT$CHL,ncol = 1),'CAR'=matrix(Input_PROSPECT$CAR,ncol = 1),'ANT'=matrix(Input_PROSPECT$ANT,ncol = 1),
                               'BROWN'=matrix(Input_PROSPECT$BROWN,ncol = 1),'EWT'=matrix(Input_PROSPECT$EWT,ncol = 1),'LMA'=matrix(Input_PROSPECT$LMA,ncol = 1),
                               'PROT'=matrix(Input_PROSPECT$PROT,ncol = 1),'CBC'=matrix(Input_PROSPECT$CBC,ncol = 1),'N'=matrix(Input_PROSPECT$N,ncol = 1),
                               'alpha'=matrix(Input_PROSPECT$alpha,ncol = 1))
  return(Input_PROSPECT)
}
