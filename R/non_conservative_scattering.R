#' computes non conservative scattering conditions
#' @param m numeric.
#' @param lai numeric. Leaf Area Index
#' @param att numeric.
#' @param sigb numeric.
#' @param ks numeric.
#' @param ko numeric.
#' @param sf numeric.
#' @param sb numeric.
#' @param vf numeric.
#' @param vb numeric.
#' @param tss numeric.
#' @param too numeric.
#'
#' @return list. tdd, rdd, tsd, rsd, tdo, rdo, rsod
#'
#' @export
non_conservative_scattering <- function(m, lai, att, sigb, ks, ko, sf, sb, vf,
                                        vb, tss, too){

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

  ps <- (sf+sb*rinf)*J1ks
  qs <- (sf*rinf+sb)*J2ks
  pv <- (vf+vb*rinf)*J1ko
  qv <- (vf*rinf+vb)*J2ko

  tdd <- (1.-rinf2)*e1/denom
  rdd <- rinf*(1.-e2)/denom
  tsd <- (ps-re*qs)/denom
  rsd <- (qs-re*ps)/denom
  tdo <- (pv-re*qv)/denom
  rdo <- (qv-re*pv)/denom

  z <- jfunc2(ks,ko,lai)
  g1 <- (z-J1ks*too)/(ko+m)
  g2 <- (z-J1ko*tss)/(ks+m)

  tv1 <- (vf*rinf+vb)*g1
  tv2 <- (vf+vb*rinf)*g2

  t1 <- tv1*(sf+sb*rinf)
  t2 <- tv2*(sf*rinf+sb)
  t3 <- (rdo*qs+tdo*ps)*rinf

  # Multiple scattering contribution to bidirectional canopy reflectance
  rsod <- (t1+t2-t3)/(1.-rinf2)
  my_list <- list("tdd" = tdd, "rdd" = rdd, "tsd" = tsd,
                  "rsd" = rsd, "tdo" = tdo, "rdo" = rdo,
                  "rsod" = rsod)
  return(my_list)
}

#' @rdname prosail-deprecated
#' @export
NonConservativeScattering <- function(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too){
  .Deprecated("non_conservative_scattering")
  non_conservative_scattering(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too)
}

