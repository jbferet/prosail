# ============================================================================== =
# prosail
# Lib_SOILSPECT.R
# ============================================================================== =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================== =
# This Library includes the SOILSPECT model
# reference: Jacquemoud, S., Baret, F., Hanocq, J.F., 1992.
# Modeling spectral and bidirectional soil reflectance. Remote Sens. Environ.
# 41, 123–132. https://doi.org/10.1016/0034-4257(92)90072-R
# ============================================================================== =

#' Computation of soil optical properties with SOILSPECT
#' (derived from Hapke) for a given single scattering albedo spectrum:

#' @param albss         = single scattering albedo
#' @param rugo          = rugosity parameter
#' @param thetav        = view zenith angle (radians)
#' @param thetas        = solar azimut angle (radians)
#' @param phi           = relative azimut angle (radians)
#' @param phase1        = phase function parameter 1
#' @param phase2        = phase function parameter 2
#' @param phase3        = phase function parameter 3
#' @param phase4        = phase function parameter 4
#' @return  bidir  = soil bidirectional reflectance
#'          dirhem = soil directional hemispherical reflectance
#'          hemhem = soil bihemispherical reflectance
#' @export

soilspect <- function(albss,rugo,thetav,thetas,phi,phase1,phase2,phase3=NULL,phase4=NULL){


  p1 <- 2.9430946*10^4
  p2 <- -1.5609100*10^5
  p3 <- 3.5659242*10^5
  p4 <- -4.5912870*10^5
  p5 <- 3.6606991*10^5
  p6 <- -1.8705815*10^5
  p7 <- 6.1339913*10^4
  p8 <- -1.2596533*10^4
  p9 <- 1.5444151*10^3
  p10 <- -1.0823175*10^2
  p11 <- 5.7691497

  if (thetav < 0){
    thetav <- -thetav
    phi <- phi +180
  }

  if (phi>=360){
    phi <- phi-360
  }

  if (phi>180){
    phi <- 360-phi
  }

  raddeg <- pi/180.
  ci <- cos(thetas*raddeg)
  ce <- cos(thetav*raddeg)
  si <- sin(thetas*raddeg)
  se <- sin(thetav*raddeg)

  #c   Computation of the Phase function : pg
  rb  <-  phase1
  rc  <-  phase2
  g1 <- ci*ce+si*se*cos(phi*raddeg)
  pg <- 1+rb*g1+(rc/2.)*(3.*g1^2.-1)

  #on ajoute deux autres param?tres
  if (!is.null(phase3) & !is.null(phase4)){
    rbprim <- phase3
    rcprim <- phase4
    g1prim <- ci*ce-si*se*cos(phi*raddeg)
    pg <- pg+rbprim*g1prim+(rcprim/2.)*(3.*g1prim^2.-1)
  }

  # Computation of the Backscattering function : bg
  bg <- 1./(1.+tan(acos(g1)/2.)/rugo)

  # Computation of the Chandrasaekar function for incidence and viewing : hi he
  g   <-  sqrt(1-albss)
  gg  <-  2.*g
  hi  <-  (1+2.*ci)/(1+gg*ci)
  he  <-  (1+2.*ce)/(1+gg*ce)

  # Computation of the bidirectional reflectance : bidir
  bidir  <-  albss/(4.*(ci+ce))*((1+bg)*pg+hi*he-1)

  # Computation of the directional hemispherical reflectance : dirhem
  dg1  <-  log(gg+1.)
  dg2  <-  log(abs(gg-1.))

  r1 <- 2.*(g+1)*ci^3*log((1.+ci)/ci)/(2.*g*ci-1)
  r2 <- (1+2.*ci)*dg1/(4.*g^2*(2.*g*ci-1.))
  r3 <- (2.*ci*(1.+g)+1.)/(2.*g)
  r4 <- ci*(4.*rb*ci^2-rc*(3.*ci^2-1.)^2-4.)*log((1.+ci)/ci)/4.
  r5 <- rb*ci*(2.*ci-1.)/2.
  r6 <- 3./8.*rc*ci*(6.*ci^3-3*ci^2-2.*ci+1.)

  dirhem  <-  (r1-r2-r3)*albss*(g-1.)/(2.*g*ci+1.)+(r4-r5+r6+1.)*albss/2

  # Computation of the bihemispherical reflectance : hemhem
  r11 <- dg1*(-p1+p2*gg-p3*gg^2+p4*gg^3-p5*gg^4+p6*gg^5-p7*gg^6
              +p8*gg^7-p9*gg^8+p10*gg^9-p11*gg^10)
  r12 <- dg2*( p1+p2*gg+p3*gg^2+p4*gg^3+p5*gg^4+p6*gg^5+p7*gg^6
               +p8*gg^7+p9*gg^8+p10*gg^9+p11*gg^10)
  r13 <- (1./13.*p1+1./12.*p2+1./11.*p3+1./10.*p4+1./9.*p5+1./8.*p6
          +1./7.*p7+1./6.*p8+1./5.*p9+1./4.*p10+1./3.*p11)*2.*gg^13
  r14 <- (1./11.*p1+1./10.*p2+1./9.*p3+1./8.*p4+1./7.*p5+1./6.*p6
          +1./5.*p7+1./4.*p8+1./3.*p9+1./2.*p10+p11)*2.*gg^11
  r15 <- (1./9.*p1+1./8.*p2+1./7.*p3+1./6.*p4+1./5.*p5+1./4.*p6
          +1./3.*p7+1./2.*p8+p9)*2.*gg^9
  r16 <- (1./7.*p1+1./6.*p2+1./5.*p3+1./4.*p4+1./3.*p5+1./2.*p6+p7)*2.*gg^7
  r17 <- (1./5.*p1+1./4.*p2+1./3.*p3+1./2.*p4+p5)*2.*gg^5
  r18 <- (1./3.*p1+1./2.*p2+p3)*2.*gg^3
  r19 <- p1*2.*gg

  rr1 <- (r11+r12+r13+r14+r15+r16+r17+r18+r19)*albss*(g^2-1.)/gg^15
  rr2 <- albss*(1.-g)*(dg1+gg*(g^2+g-1.))/(g*gg^3)
  rr3 <- albss*(1.-g)*dg1*(4.*g+(g-1.)*dg1+(g+1.)*dg2)/gg^5

  rr41 <- -4.*(1./13.*p1+1./12.*p2+1./11.*p3+1./10.*p4+1./9.*p5
               +1./8.*p6+1./7.*p7+1./6.*p8+1./5.*p9+1./4.*p10+1./3.*p11)
  rr42 <- 4.*rb*(1./15.*p1+1./14.*p2+1./13.*p3+1./12.*p4+1./11.*p5
                 +1./10.*p6+1./9.*p7+1./8.*p8+1./7.*p9+1./6.*p10+1./5.*p11)
  rr43 <- -rc*(228./1105.*p1+73./336.*p2+164./715.*p3+17./70.*p4
               +332./1287.*p5+11./40.*p6+68./231.*p7+19./60.*p8+12./35.*p9
               +3./8.*p10+44./105.*p11)

  rr4 <- (rr41+rr42+rr43)*albss/8.
  rr5 <- -albss*rb/24.
  rr6 <- 7*albss*rc/160.
  rr7 <- albss/4.

  hemhem <- 2.*(rr1+rr2+rr3+rr4+rr5+rr6+rr7)
  mylist <- list('bidir'=bidir,'dirhem'=dirhem,'hemhem'=hemhem)
  return(mylist)
}
