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
  SZA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SZA <- c(SZA,as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,
                                     split = ' ')[[1]]))
  }
  SZA <- stats::na.omit(SZA)

  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SAA <- c(SAA,as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,
                                     split = ' ')[[1]]))
  }
  SAA <- stats::na.omit(SAA)

  # VZA
  VZA <- c()
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    values <- band_detector[i]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
    for (j in seq_len(length(values))){
      VZA <- c(VZA,as.numeric(strsplit(x = values[j]$VALUES,
                                       split = ' ')[[1]]))
    }
  }
  VZA <- stats::na.omit(VZA)

  # VAA
  VAA <- c()
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    values <- band_detector[i]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
    for (j in seq_len(length(values))){
      VAA <- c(VAA,as.numeric(strsplit(x = values[j]$VALUES,
                                       split = ' ')[[1]]))
    }
  }
  VAA <- stats::na.omit(VAA)
  return(list('SAA' = SAA, 'SZA' = SZA, 'VAA' = VAA, 'VZA' = VZA))
}
