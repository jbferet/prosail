#' computes conservative scattering conditions
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
conservative_scattering <- function(m, lai, att, sigb, ks, ko, sf, sb, vf, vb,
                                    tss, too){

  # Near or complete conservative scattering
  J4 <- Jfunc4(m,lai)
  amsig <- att-sigb
  apsig <- att+sigb
  rtp <- (1-amsig*J4)/(1+amsig*J4)
  rtm <- (-1+apsig*J4)/(1+apsig*J4)
  rdd <- 0.5*(rtp+rtm)
  tdd <- 0.5*(rtp-rtm)

  dns <- ks*ks-m*m
  dno <- ko*ko-m*m
  cks <- (sb*(ks-att) - sf*sigb)/dns
  cko <- (vb*(ko-att) - vf*sigb)/dno
  dks <- (-sf*(ks+att) - sb*sigb)/dns
  dko <- (-vf*(ko+att) - vb*sigb)/dno
  ho <- (sf*cko+sb*dko)/(ko+ks)

  rsd <- cks*(1-tss*tdd) - dks*rdd
  rdo <- cko*(1-too*tdd) - dko*rdd
  tsd <- dks*(tss-tdd) - cks*tss*rdd
  tdo <- dko*(too-tdd) - cko*too*rdd
  # Multiple scattering contribution to bidirectional canopy reflectance
  rsod <- ho*(1-tss*too) - cko*tsd*too - dko*rsd

  my_list <- list("tdd" = tdd, "rdd" = rdd, "tsd" = tsd, "rsd" = rsd,
                  "tdo" = tdo, "rdo" = rdo, "rsod" = rsod)
  return(my_list)
}
