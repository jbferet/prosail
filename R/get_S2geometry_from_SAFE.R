#' This function returns geometry of acquisition for S2 image processed with
#' Sen2Cor
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#' @param s2xml list. produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @importFrom stats na.omit
#' @export

get_S2geometry_from_SAFE <- function(s2xml){

  Distrib_SunAngle <- list()
  # SZA
  Distrib_SunAngle$Zenith <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Zenith$Values_List
  SZA <- list()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SZA[[i]] <- as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,
                                    split = ' ')[[1]])
  }
  SZA <- do.call(what = rbind, args = SZA)
  # SZA <- stats::na.omit(SZA)

  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- list()
  for (i in seq_len(length(Distrib_SunAngle$Azimuth))){
    SAA[[i]] <- as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,
                                    split = ' ')[[1]])
  }
  SAA <- do.call(what = rbind, args = SAA)

  # VZA
  VZA <- list()
  ii <- 0
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    ii <- ii + 1
    vza_d <- band_detector[i]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
    for (j in seq_len(length(Distrib_SunAngle$Zenith))){
      vza_d[[j]] <- as.numeric(strsplit(x = vza_d[[j]],
                                      split = ' ')[[1]])
    }
    VZA[[ii]] <- do.call(what = rbind, args = vza_d)
  }
  # get mean view zenith angle for all detectors
  VZA <- apply(simplify2array(VZA), 1:2, mean, na.rm=T)

  # VAA
  VAA <- list()
  ii <- 0
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    ii <- ii + 1
    vaa_d <- band_detector[i]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
    for (j in seq_len(length(Distrib_SunAngle$Azimuth))){
      vaa_d[[j]] <- as.numeric(strsplit(x = vaa_d[[j]],
                                        split = ' ')[[1]])
    }
    VAA[[ii]] <- do.call(what = rbind, args = vaa_d)
  }
  # get mean view azimuth angle for all detectors
  VAA <- apply(simplify2array(VAA), 1:2, mean, na.rm=T)
  return(list('SAA' = SAA, 'SZA' = SZA, 'VAA' = VAA, 'VZA' = VZA))
}
