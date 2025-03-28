#' returns geometry of acquisition for S2 image processed with MAJA
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#' @param s2xml list. produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @importFrom stats na.omit
#' @export

get_s2_geometry_from_THEIA <- function(s2xml){

  Distrib_SunAngle <- list()
  # SZA
  Distrib_SunAngle$Zenith <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Zenith$Values_List
  SZA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith)))
    SZA <- c(SZA, as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,
                                      split = ' ')[[1]]))
  # SZA <- stats::na.omit(SZA)

  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith)))
    SAA <- c(SAA,as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,
                                     split = ' ')[[1]]))
  # SAA <- stats::na.omit(SAA)

  # VZA
  VZA <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
      for (k in seq_len(length(values)))
        VZA <- c(VZA,as.numeric(strsplit(x = values[k]$VALUES,
                                         split = ' ')[[1]]))
    }
  }
  # VZA <- stats::na.omit(VZA)

  # VAA
  VAA <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
      for (k in seq_len(length(values)))
        VAA <- c(VAA,as.numeric(strsplit(x = values[k]$VALUES,
                                         split = ' ')[[1]]))
    }
  }
  # VAA <- stats::na.omit(VAA)
  return(list('SAA' = SAA, 'SZA' = SZA, 'VAA' = VAA, 'VZA' = VZA))
}
