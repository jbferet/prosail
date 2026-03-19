#' 4SAIL2 simulation based on a set of combinations of input parameters
#' @param leaf_green dataframe. leaf optical properties #1 (e.g. green veg)
#' @param leaf_brown dataframe. leaf optical properties #2 (e.g. brown veg)
#' @param type_lidf numeric. Type of leaf inclination distribution function
#' @param lidf_a numeric.
#' if type_lidf ==1, controls the average leaf slope
#' if type_lidf ==2, controls the average leaf angle
#' @param lidf_b numeric.
#' if type_lidf ==1, unused
#' if type_lidf ==2, controls the distribution's bimodality
#' @param lai numeric. Leaf Area Index
#' @param hotspot numeric. Hot Spot parameter = ratio of the correlation length of
#' leaf projections in the horizontal plane and the canopy height
#' (doi:10.1016/j.rse.2006.12.013)
#' @param tts numeric. Sun zeith angle
#' @param tto numeric. Observer zeith angle
#' @param psi numeric. Azimuth Sun / Observer
#' @param rsoil numeric. Soil reflectance
#' @param fraction_brown numeric. Fraction of brown leaf area
#' @param diss numeric. Layer dissociation factor
#' @param cv numeric. vertical crown cover percentage
#' = % ground area covered with crowns as seen from nadir direction
#' @param zeta numeric. Tree shape factor
#' = ratio of crown diameter to crown height
#'
#' @return list. rdot,rsot,rddt,rsdt
#' rdot: hemispherical-directional reflectance factor in viewing direction
#' rsot: bi-directional reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rddt: bi-hemispherical reflectance factor
#' fcover: Fraction of green Vegetation Cover (= 1 - beam transmittance in the
#' target-view path)
#' abs_dir: canopy absorptance for direct solar incident flux
#' abs_hem: canopy absorptance for hemispherical diffuse incident flux
#' rsdstar: contribution of direct solar incident flux to albedo
#' rddstar: contribution of hemispherical diffuse incident flux to albedo

#' @export

fourSAIL2  <- function(leaf_green, leaf_brown, type_lidf = 2, lidf_a = 60,
                       lidf_b = NULL, lai = 3, hotspot = 0.1, tts = 30, tto = 0,
                       psi = 60, rsoil = NULL, fraction_brown = 0.5, diss = 0.5,
                       cv = 1, zeta = 1){

  #	This version does not include non-Lambertian soil properties.
  #	original codes do, and only need to add the following variables as input
  rddsoil <- rdosoil <- rsdsoil <- rsosoil <- rsoil

  #	Geometric quantities
  rd <- pi/180

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

  if (lai<0){
    message('Please define positive LAI value')
    rddt <- rsdt <- rdot <- rsost <- rsot <- rsoil
    alfast <- alfadt <- 0*rsoil
  } else if (lai==0){
    tss <- too <- tsstoo <- tdd <- 1.0
    rdd <- rsd <- tsd <- rdo <- tdo <- 0.0
    rso <- rsos <- rsod <- rsodt <- 0.0
    rddt <- rsdt <- rdot <- rsost <- rsot <- rsoil
    alfast <- alfadt <- 0*rsoil
  } else if (lai>0){
    cts <- cos(rd*tts)
    cto <- cos(rd*tto)
    ctscto <- cts*cto
    tants <- tan(rd*tts)
    tanto <- tan(rd*tto)
    cospsi <- cos(rd*psi)
    dso <- sqrt(tants*tants+tanto*tanto-2.0*tants*tanto*cospsi)

    # Clumping effects
    Cs <- Co <- 1.0
    if (cv<=1.0){
      Cs <- 1.0-(1.0-cv)^(1.0/cts)
      Co <- 1.0-(1.0-cv)^(1.0/cto)
    }
    overlap <- 0.0
    if (zeta>0.0)
      overlap <- min(Cs*(1.0-Co),Co*(1.0-Cs))*exp(-dso/zeta)
    Fcd <- Cs*Co+overlap
    Fcs <- (1.0-Cs)*Co-overlap
    Fod <- Cs*(1.0-Co)-overlap
    Fos <- (1.0-Cs)*(1.0-Co)+overlap
    Fcdc <- 1.0-(1.0-Fcd)^(0.5/cts+0.5/cto)

    #	Part depending on diss, fraction_brown, and leaf optical properties
    #	First save the input fraction_brown as the old fraction_brown, as the
    # following change is only artificial
    # Better define an fraction_brown that is actually used: fb, so that the
    # input is not modified!

    fb <- fraction_brown
    # if only green leaves
    if (fraction_brown==0.0){
      fb <- 0.5
      leaf_brown$Reflectance <- leaf_green$Reflectance
      leaf_brown$Transmittance <- leaf_green$Transmittance
    }
    if (fraction_brown==1.0){
      fb <- 0.5
      leaf_green$Reflectance <- leaf_brown$Reflectance
      leaf_green$Transmittance <- leaf_brown$Transmittance
    }
    s <- (1.0-diss)*fb*(1.0-fb)
    # rho1 & tau1 : green foliage
    # rho2 & tau2 : brown foliage (bottom layer)
    rho1 <- ((1-fb-s)*leaf_green$Reflectance+s*leaf_brown$Reflectance)/(1-fb)
    tau1 <- ((1-fb-s)*leaf_green$Transmittance+s*leaf_brown$Transmittance)/(1-fb)
    rho2 <- (s*leaf_green$Reflectance+(fb-s)*leaf_brown$Reflectance)/fb
    tau2 <- (s*leaf_green$Transmittance+(fb-s)*leaf_brown$Transmittance)/fb

    # angular distance, compensation of shadow length
    #	Calculate geometric factors associated with extinction and scattering
    #	Initialise sums
    ks <- ko <- bf <- sob <- sof <- 0

    # Weighted sums over LIDF

    for (i in seq_len(length(litab))){
      ttl <- litab[i]
      ctl <- cos(rd*ttl)
      # SAIL volscatt function gives interception coefficients
      # and two portions of the volume scattering phase function to be
      # multiplied by rho and tau, respectively
      res_volscatt <- volscatt(tts,tto,psi,ttl)
      chi_s <- res_volscatt$chi_s
      chi_o <- res_volscatt$chi_o
      frho <- res_volscatt$frho
      ftau <- res_volscatt$ftau
      # Extinction coefficients
      ksli <- chi_s/cts
      koli <- chi_o/cto
      # Area scattering coefficient fractions
      sobli <- frho*pi/ctscto
      sofli <- ftau*pi/ctscto
      bfli <- ctl*ctl
      ks <- ks+ksli*lidf[i]
      ko <- ko+koli*lidf[i]
      bf <- bf+bfli*lidf[i]
      sob <- sob+sobli*lidf[i]
      sof <- sof+sofli*lidf[i]
    }
    # Geometric factors to be used later in combination with rho and tau
    sdb <- 0.5*(ks+bf)
    sdf <- 0.5*(ks-bf)
    dob <- 0.5*(ko+bf)
    dof <- 0.5*(ko-bf)
    ddb <- 0.5*(1.+bf)
    ddf <- 0.5*(1.-bf)

    # LAIs in two layers
    lai1 <- (1-fb)*lai
    lai2 <- fb*lai

    tss <- exp(-ks*lai)
    ck <- exp(-ks*lai1)
    alf <- 1e6
    if (hotspot > 0)
      alf <- (dso/hotspot)*2./(ks+ko)
    if (alf>200.0)
      alf<- 200.0     # inserted H. Bach 1/3/04
    if (alf==0.0){
      # The pure hotspot
      tsstoo <- tss
      s1 <- (1-ck)/(ks*lai)
      s2 <- (ck-tss)/(ks*lai)
    } else {
      # Outside the hotspot
      fhot <- lai*sqrt(ko*ks)
      # Integrate 2 layers by exponential simpson method in 20 steps
      # the steps are arranged according to equal partitioning
      # of the derivative of the joint probability function
      x1 <- y1 <- 0.0
      f1 <- 1.0
      ca <- exp(alf*(fb-1.0))
      fint <- (1.0-ca)*.05
      s1 <- 0.0
      for (istep in 1:20){
        if (istep<20){
          x2 <- -log(1.-istep*fint)/alf
        } else {
          x2 <- 1.-fb
        }
        y2 <- -(ko+ks)*lai*x2+fhot*(1.0-exp(-alf*x2))/alf
        f2 <- exp(y2)
        s1 <- s1+(f2-f1)*(x2-x1)/(y2-y1)
        x1 <- x2
        y1 <- y2
        f1 <- f2
      }
      fint <- (ca-exp(-alf))*.05
      s2 <- 0.0
      for (istep in 1:20){
        if (istep<20){
          x2 <- -log(ca-istep*fint)/alf
        } else {
          x2 <- 1.0
        }
        y2 <- -(ko+ks)*lai*x2+fhot*(1.0-exp(-alf*x2))/alf
        f2 <- exp(y2)
        s2 <- s2+(f2-f1)*(x2-x1)/(y2-y1)
        x1 <- x2
        y1 <- y2
        f1 <- f2
      }
      tsstoo <- f1
    }

    # Calculate reflectances and transmittances
    # Bottom layer
    tss <- exp(-ks*lai2)
    too <- exp(-ko*lai2)
    sb <- sdb*rho2+sdf*tau2
    sf <- sdf*rho2+sdb*tau2

    vb <- dob*rho2+dof*tau2
    vf <- dof*rho2+dob*tau2

    w2 <- sob*rho2+sof*tau2

    sigb <- ddb*rho2+ddf*tau2
    sigf <- ddf*rho2+ddb*tau2
    att <- 1.0-sigf
    m2 <- (att+sigb)*(att-sigb)
    m2[m2<0] <- 0
    m <- sqrt(m2)
    which_ncs <- which(m>0.01)
    which_cs <- which(m<=0.01)

    tdd <- rdd <- tsd <- rsd <- tdo <- rdo <- 0*m
    rsod <- 0*m
    if (length(which_ncs)>0){
      res_ncs <- non_conservative_scattering(m[which_ncs], lai2, att[which_ncs],
                                             sigb[which_ncs], ks, ko,
                                             sf[which_ncs], sb[which_ncs],
                                             vf[which_ncs], vb[which_ncs],
                                             tss, too)
      tdd[which_ncs] <- res_ncs$tdd
      rdd[which_ncs] <- res_ncs$rdd
      tsd[which_ncs] <- res_ncs$tsd
      rsd[which_ncs] <- res_ncs$rsd
      tdo[which_ncs] <- res_ncs$tdo
      rdo[which_ncs] <- res_ncs$rdo
      rsod[which_ncs] <- res_ncs$rsod
    }
    if (length(which_cs)>0){
      res_cs <- conservative_scattering(m[which_cs], lai2, att[which_cs],
                                        sigb[which_cs], ks, ko, sf[which_cs],
                                        sb[which_cs], vf[which_cs],
                                        vb[which_cs], tss, too)
      tdd[which_cs] <- res_cs$tdd
      rdd[which_cs] <- res_cs$rdd
      tsd[which_cs] <- res_cs$tsd
      rsd[which_cs] <- res_cs$rsd
      tdo[which_cs] <- res_cs$tdo
      rdo[which_cs] <- res_cs$rdo
      rsod[which_cs] <- res_cs$rsod
    }

    # Set background properties equal to those of the bottom layer on black soil
    rddb <- rdd
    rsdb <- rsd
    rdob <- rdo
    rsodb <- rsod
    tddb <- tdd
    tsdb <- tsd
    tdob <- tdo
    toob <- too
    tssb <- tss
    # Top layer
    tss <- exp(-ks*lai1)
    too <- exp(-ko*lai1)

    sb <- sdb*rho1+sdf*tau1
    sf <- sdf*rho1+sdb*tau1

    vb <- dob*rho1+dof*tau1
    vf <- dof*rho1+dob*tau1

    w1 <- sob*rho1+sof*tau1

    sigb <- ddb*rho1+ddf*tau1
    sigf <- ddf*rho1+ddb*tau1
    att <- 1.0-sigf

    m2 <- (att+sigb)*(att-sigb)
    m2[m2<0] <- 0
    m <- sqrt(m2)
    which_ncs <- which(m>0.01)
    which_cs <- which(m<=0.01)

    tdd <- rdd <- tsd <- rsd <- tdo <- rdo <- 0*m
    rsod <- 0*m
    if (length(which_ncs)>0){
      res_ncs <- non_conservative_scattering(m[which_ncs], lai1, att[which_ncs],
                                             sigb[which_ncs], ks, ko,
                                             sf[which_ncs], sb[which_ncs],
                                             vf[which_ncs], vb[which_ncs],
                                             tss, too)
      tdd[which_ncs] <- res_ncs$tdd
      rdd[which_ncs] <- res_ncs$rdd
      tsd[which_ncs] <- res_ncs$tsd
      rsd[which_ncs] <- res_ncs$rsd
      tdo[which_ncs] <- res_ncs$tdo
      rdo[which_ncs] <- res_ncs$rdo
      rsod[which_ncs] <- res_ncs$rsod
    }
    if (length(which_cs)>0){
      res_cs <- conservative_scattering(m[which_cs], lai1, att[which_cs],
                                        sigb[which_cs], ks, ko, sf[which_cs],
                                        sb[which_cs], vf[which_cs],
                                        vb[which_cs], tss, too)
      tdd[which_cs] <- res_cs$tdd
      rdd[which_cs] <- res_cs$rdd
      tsd[which_cs] <- res_cs$tsd
      rsd[which_cs] <- res_cs$rsd
      tdo[which_cs] <- res_cs$tdo
      rdo[which_cs] <- res_cs$rdo
      rsod[which_cs] <- res_cs$rsod
    }

    # Combine with bottom layer reflectances and transmittances (adding method)
    rn <- 1.0-rdd*rddb
    tup <- (tss*rsdb+tsd*rddb)/rn
    tdn <- (tsd+tss*rsdb*rdd)/rn
    rsdt <- rsd+tup*tdd
    rdot <- rdo+tdd*(rddb*tdo+rdob*too)/rn
    rsodt <- rsod+(tss*rsodb+tdn*rdob)*too+tup*tdo

    rsost <- (w1*s1+w2*s2)*lai

    rsot <- rsost+rsodt

    # Diffuse reflectances at the top and the bottom are now different
    rddt_t <- rdd+tdd*rddb*tdd/rn
    rddt_b <- rddb+tddb*rdd*tddb/rn

    # Transmittances of the combined canopy layers
    tsst <- tss*tssb
    toot <- too*toob
    tsdt <- tss*tsdb+tdn*tddb
    tdot <- tdob*too+tddb*(tdo+rdd*rdob*too)/rn
    tddt <- tdd*tddb/rn

    # Apply clumping effects to vegetation layer
    rddcb <- cv*rddt_b
    rddct <- cv*rddt_t
    tddc <- 1-cv+cv*tddt
    rsdc <- Cs*rsdt
    tsdc <- Cs*tsdt
    rdoc <- Co*rdot
    tdoc <- Co*tdot
    tssc <- 1-Cs+Cs*tsst
    tooc <- 1-Co+Co*toot

    # New weight function Fcdc for crown contribution (W. Verhoef, 22-05-08)
    rsoc <- Fcdc*rsot
    tssooc <- Fcd*tsstoo+Fcs*toot+Fod*tsst+Fos

    # Canopy absorptance for black background (W. Verhoef, 02-03-04)
    alfas <- 1.-tssc-tsdc-rsdc
    alfad <- 1.-tddc-rddct

    # Add the soil background
    rn <- 1-rddcb*rddsoil
    tup <- (tssc*rsdsoil+tsdc*rddsoil)/rn
    tdn <- (tsdc+tssc*rsdsoil*rddcb)/rn

    rddt <- rddct+tddc*rddsoil*tddc/rn
    rsdt <- rsdc+tup*tddc
    rdot <- rdoc+tddc*(rddsoil*tdoc+rdosoil*tooc)/rn
    rsot <- rsoc+tssooc*rsosoil+tdn*rdosoil*tooc+tup*tdoc

    # Effect of soil background on canopy absorptances (W. Verhoef, 02-03-04)
    alfast <- alfas+tup*alfad
    alfadt <- alfad*(1.+tddc*rddsoil/rn)

    # # compute directional and hemispherical absorbances
    # abs_dir <- 1 - rsdt - ((1-rsoil)*tss) - (1-rsoil)*((tss*rsoil*rdd)+tsd)/dn
    # abs_hem <- 1 - rddt - ((1-rsoil)*tdd) - (1-rsoil)*(tdd*rdd*rsoil)/dn

    # # compute absorbances of the isolated canopy (from J. Gomez Dans)
    # alfas <- 1.0 - tss - tsd - rsd  # direct flux
    # alfad <- 1.0 - tdd - rdd  # diffuse
    # alfasx <- alfas + (rsoil * (tss + tsd) / dn) * alfad
    # alfadx <- alfad + ((tdd * rsoil) / dn) * alfas

    # compute Albedo (from J. Gomez Dans)
    rsdstar <- rsd + (tss + tsd) * rsoil * tdd / rn
    rddstar <- rdd + (tdd * tdd * rsoil) / rn
  }
  my_list <- list('rdot' = rdot, 'rsot' = rsot, 'rddt' = rddt, 'rsdt' =rsdt,
                  'fcover' = 1 - too, 'abs_dir' = alfast, 'abs_hem' = alfadt,
                  'rsdstar' = rsdstar, 'rddstar' = rddstar)
  return(my_list)
}
