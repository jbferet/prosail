#' Performs PROSAIL simulation based on a set of combinations of
#' input parameters
#' @param lop list. leaf optical properties
#' (reflectance and transmittance) and corresponding spectral bands
#' @param type_lidf numeric. Type of leaf inclination distribution function
#' @param lidf_a numeric.
#' if type_lidf ==1, controls the average leaf slope
#' if type_lidf ==2, controls the average leaf angle
#' @param lidf_b numeric.
#' if type_lidf ==1, controls the distribution's bimodality
#' if type_lidf ==2, unused
#' @param lai numeric. Leaf Area Index
#' @param hotspot numeric. Hot Spot parameter
#' @param tts numeric. Sun zeith angle
#' @param tto numeric. Observer zeith angle
#' @param psi numeric. Azimuth Sun / Observer
#' @param rsoil numeric. Soil reflectance
#'
#' @return list. rdot,rsot,rddt,rsdt
#' rdot: hemispherical-directional reflectance factor in viewing direction
#' rsot: bi-directional reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rddt: bi-hemispherical reflectance factor
#' fcover: Fraction of green Vegetation Cover (= 1 - beam transmittance in
#' the target-view path)
#' abs_dir: canopy absorptance for direct solar incident flux
#' abs_hem: canopy absorptance for hemispherical diffuse incident flux
#' rsdstar: contribution of direct solar incident flux to albedo
#' rddstar: contribution of hemispherical diffuse incident flux to albedo
#' @export

fourSAIL  <- function(lop, type_lidf = 2, lidf_a = 60, lidf_b = NULL,
                      lai = 3, hotspot = 0.1, tts = 30, tto = 0, psi = 60,
                      rsoil = NULL){

  ############################ #
  #	LEAF OPTICAL PROPERTIES	##
  ############################ #
  rho <- lop$reflectance
  tau <- lop$transmittance

  #	Geometric quantities
  rd <- pi/180
  cts <- cos(rd*tts)
  cto <- cos(rd*tto)
  ctscto <- cts*cto
  tants <- tan(rd*tts)
  tanto <- tan(rd*tto)
  cospsi <- cos(rd*psi)
  dso <- sqrt(tants*tants+tanto*tanto-2.*tants*tanto*cospsi)

  #	Generate leaf angle distribution from average leaf angle (ellipsoidal)
  # or (a,b) parameters
  if (type_lidf==1){
    foliar_distrib <- dladgen(lidf_a,lidf_b)
    lidf <- foliar_distrib$lidf
    litab <- foliar_distrib$litab

  } else if (type_lidf==2){
    foliar_distrib <- campbell(lidf_a)
    lidf <- foliar_distrib$lidf
    litab <- foliar_distrib$litab
  }

  # angular distance, compensation of shadow length
  #	Calculate geometric factors associated with extinction and scattering
  #	Initialise sums
  ks <- 0
  ko <- 0
  bf <- 0
  sob <- 0
  sof <- 0

  #	Weighted sums over LIDF
  na <- length(litab)
  for (i in 1:na){
    ttl <- litab[i]	    # leaf inclination discrete values
    ctl <- cos(rd*ttl)
    #	SAIL volume scattering phase function gives interception and portions to
    #	be multiplied by rho and tau
    res_volscatt <- volscatt(tts,tto,psi,ttl)
    chi_s <- res_volscatt$chi_s
    chi_o <- res_volscatt$chi_o
    frho <- res_volscatt$frho
    ftau <- res_volscatt$ftau

    # **************************************************************************
    #                   SUITS SYSTEM COEFFICIENTS
    #
    # ks : Extinction coefficient for direct solar flux
    # ko : Extinction coefficient for direct observed flux
    # att : Attenuation coefficient for diffuse flux
    # sigb : Backscattering coefficient of the diffuse downward flux
    # sigf : Forwardscattering coefficient of the diffuse upward flux
    # sf : Scattering coef of direct solar flux for downward diffuse flux
    # sb : Scattering coef of  direct solar flux for upward diffuse flux
    # vf : Scattering coef of upward diffuse flux in the observed direction
    # vb : Scattering coef of downward diffuse flux in the observed direction
    # w : Bidirectional scattering coefficient
    # **************************************************************************

    #	Extinction coefficients
    ksli <- chi_s/cts
    koli <- chi_o/cto

    #	Area scattering coefficient fractions
    sobli <- frho*pi/ctscto
    sofli <- ftau*pi/ctscto
    bfli <- ctl*ctl
    ks <- ks+ksli*lidf[i]
    ko <- ko+koli*lidf[i]
    bf <- bf+bfli*lidf[i]
    sob <- sob+sobli*lidf[i]
    sof <- sof+sofli*lidf[i]
  }

  #	Geometric factors to be used later with rho and tau
  sdb <- 0.5*(ks+bf)
  sdf <- 0.5*(ks-bf)
  dob <- 0.5*(ko+bf)
  dof <- 0.5*(ko-bf)
  ddb <- 0.5*(1.+bf)
  ddf <- 0.5*(1.-bf)

  #	Here rho and tau come in
  sigb <- ddb*rho+ddf*tau
  sigf <- ddf*rho+ddb*tau
  att <- 1-sigf
  m2 <- (att+sigb)*(att-sigb)
  m2[which(m2<=0)] <- 0
  m <- sqrt(m2)

  sb <- sdb*rho+sdf*tau
  sf <- sdf*rho+sdb*tau
  vb <- dob*rho+dof*tau
  vf <- dof*rho+dob*tau
  w <- sob*rho+sof*tau

  #	Here the LAI comes in
  #   Outputs for the case LAI = 0
  if (lai<0){
    tss <- 1          # beam transmittance in the sun-target path.
    too <- 1          # beam transmittance in the target-view path.
    tsstoo <- 1       # beam transmittance in the sun-target-view path.
    rdd <- 0          # canopy bihemispherical reflectance factor.
    tdd <- 1          # canopy bihemispherical transmittance factor.
    rsd <- 0          # canopy directional-hemispherical reflectance factor.
    tsd <- 0          # canopy directional-hemispherical transmittance factor.
    rdo <- 0          # canopy hemispherical-directional reflectance factor.
    tdo <- 0          # canopy hemispherical-directional transmittance factor.
    rso <- 0          # canopy bidirectional reflectance factor.
    rsos <- 0         # single scattering contribution to rso.
    rsod <- 0         # multiple scattering contribution to rso.

    rddt <- rsoil     # surface bihemispherical reflectance factor.
    rsdt <- rsoil     # surface directional-hemispherical reflectance factor.
    rdot <- rsoil     # surface hemispherical-directional reflectance factor.
    rsodt <- 0*rsoil  # reflectance factor.
    rsost <- rsoil    # reflectance factor.
    rsot <- rsoil     # surface bidirectional reflectance factor.
  } else {
    #	Other cases (LAI > 0)
    e1 <- exp(-m*lai)
    e2 <- e1*e1
    rinf <- (att-m)/sigb
    rinf2 <- rinf*rinf
    re <- rinf*e1
    denom <- 1.-rinf2*e2

    J1ks <- jfunc1(ks,m,lai)
    J2ks <- jfunc2(ks,m,lai)
    J1ko <- jfunc1(ko,m,lai)
    J2ko <- jfunc2(ko,m,lai)

    Ps <- (sf+sb*rinf)*J1ks
    Qs <- (sf*rinf+sb)*J2ks
    Pv <- (vf+vb*rinf)*J1ko
    Qv <- (vf*rinf+vb)*J2ko

    rdd <- rinf*(1.-e2)/denom
    tdd <- (1.-rinf2)*e1/denom
    tsd <- (Ps-re*Qs)/denom
    rsd <- (Qs-re*Ps)/denom
    tdo <- (Pv-re*Qv)/denom
    rdo <- (Qv-re*Pv)/denom

    tss <- exp(-ks*lai)
    too <- exp(-ko*lai)
    z <- jfunc3(ks,ko,lai)
    g1 <- (z-J1ks*too)/(ko+m)
    g2 <- (z-J1ko*tss)/(ks+m)

    Tv1 <- (vf*rinf+vb)*g1
    Tv2 <- (vf+vb*rinf)*g2
    T1 <- Tv1*(sf+sb*rinf)
    T2 <- Tv2*(sf*rinf+sb)
    T3 <- (rdo*Qs+tdo*Ps)*rinf

    #	Multiple scattering contribution to bidirectional canopy reflectance
    rsod <- (T1+T2-T3)/(1.-rinf2)

    #	Treatment of the hotspot-effect
    alf <- 1e6
    #	Apply correction 2/(K+k) suggested by F.-M. Breon
    if (hotspot>0)
      alf <- (dso/hotspot)*2./(ks+ko)
    # inserted H. Bach 1/3/04
    if (alf>200)
      alf <- 200
    if (alf==0){
      #	The pure hotspot - no shadow
      tsstoo <- tss
      sumint <- (1-tss)/(ks*lai)
    } else {
      #	Outside the hotspot
      fhot <- lai*sqrt(ko*ks)
      #	Integrate by exponential Simpson method in 20 steps
      #	the steps are arranged according to equal partitioning
      #	of the slope of the joint probability function
      x1 <- y1 <- sumint <- 0
      f1 <- 1
      fint <- (1.-exp(-alf))*0.05
      for (i in 1:20){
        if (i<20){
          x2 <- -log(1.-i*fint)/alf
        } else {
          x2 <- 1
        }
        y2 <- -(ko+ks)*lai*x2+fhot*(1.-exp(-alf*x2))/alf
        f2 <- exp(y2)
        sumint <- sumint+(f2-f1)*(x2-x1)/(y2-y1)
        x1 <- x2
        y1 <- y2
        f1 <- f2
      }
      tsstoo <- f1
    }
    #	Bidirectional reflectance
    #	Single scattering contribution
    rsos <- w*lai*sumint
    #	Total canopy contribution
    rso <- rsos+rsod
    #	Interaction with the soil
    dn <- 1.-rsoil*rdd
    # rddt: bi-hemispherical reflectance factor
    rddt <- rdd+tdd*rsoil*tdd/dn
    # rsdt: directional-hemispherical reflectance factor for solar incident flux
    rsdt <- rsd+(tsd+tss)*rsoil*tdd/dn
    # rdot: hemispherical-directional reflectance factor in viewing direction
    rdot <- rdo+tdd*rsoil*(tdo+too)/dn
    # rsot: bi-directional reflectance factor
    rsodt <- rsod+((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn
    rsost <- rsos+tsstoo*rsoil
    rsot <- rsost+rsodt

    # compute directional and hemispherical absorbances
    abs_dir <- 1 - rsdt - ((1-rsoil)*tss) - (1-rsoil)*((tss*rsoil*rdd)+tsd)/dn
    abs_hem <- 1 - rddt - ((1-rsoil)*tdd) - (1-rsoil)*(tdd*rdd*rsoil)/dn

    # # compute absorbances of the isolated canopy (from J. Gomez Dans)
    # alfas <- 1.0 - tss - tsd - rsd  # direct flux
    # alfad <- 1.0 - tdd - rdd  # diffuse
    # alfasx <- alfas + (rsoil * (tss + tsd) / dn) * alfad
    # alfadx <- alfad + ((tdd * rsoil) / dn) * alfas

    # compute Albedo (from J. Gomez Dans)
    rsdstar <- rsd + (tss + tsd) * rsoil * tdd / dn
    rddstar <- rdd + (tdd * tdd * rsoil) / dn
  }
  my_list <- list('rdot' = rdot, 'rsot' = rsot, 'rddt' = rddt, 'rsdt' = rsdt,
                  'fcover' = 1 - too, 'abs_dir' = abs_dir, 'abs_hem' = abs_hem,
                  'rsdstar' = rsdstar, 'rddstar' = rddstar)
  return(my_list)
}
