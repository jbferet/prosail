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

  J1ks <- Jfunc1(ks,m,lai)
  J2ks <- Jfunc2(ks,m,lai)
  J1ko <- Jfunc1(ko,m,lai)
  J2ko <- Jfunc2(ko,m,lai)

  Ps <- (sf+sb*rinf)*J1ks
  Qs <- (sf*rinf+sb)*J2ks
  Pv <- (vf+vb*rinf)*J1ko
  Qv <- (vf*rinf+vb)*J2ko

  tdd <- (1.-rinf2)*e1/denom
  rdd <- rinf*(1.-e2)/denom
  tsd <- (Ps-re*Qs)/denom
  rsd <- (Qs-re*Ps)/denom
  tdo <- (Pv-re*Qv)/denom
  rdo <- (Qv-re*Pv)/denom

  z <- Jfunc2(ks,ko,lai)
  g1 <- (z-J1ks*too)/(ko+m)
  g2 <- (z-J1ko*tss)/(ks+m)

  Tv1 <- (vf*rinf+vb)*g1
  Tv2 <- (vf+vb*rinf)*g2

  T1 <- Tv1*(sf+sb*rinf)
  T2 <- Tv2*(sf*rinf+sb)
  T3 <- (rdo*Qs+tdo*Ps)*rinf

  # Multiple scattering contribution to bidirectional canopy reflectance
  rsod <- (T1+T2-T3)/(1.-rinf2)
  my_list <- list("tdd" = tdd, "rdd" = rdd, "tsd" = tsd,
                  "rsd" = rsd, "tdo" = tdo, "rdo" = rdo,
                  "rsod" = rsod)
  return(my_list)
}
