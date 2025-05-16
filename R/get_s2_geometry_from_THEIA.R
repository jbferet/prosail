#' returns geometry of acquisition for S2 image processed with MAJA
#' - sza = list of sun zenith angle
#' - saa = list of sun azimuth angle
#' - vza = list of viewer zenith angle
#' - vaa = list of viewer azimuth angle
#' @param s2xml list. produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (sza, saa, vza, vaa)
#' @importFrom stats na.omit
#' @export

get_s2_geometry_from_THEIA <- function(s2xml){

  distrib_sun_angle <- list()
  # sza
  distrib_sun_angle$zenith <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Zenith$Values_List
  sza <- c()
  for (i in seq_len(length(distrib_sun_angle$zenith)))
    sza <- c(sza, as.numeric(strsplit(x = distrib_sun_angle$zenith[i]$VALUES,
                                      split = ' ')[[1]]))
  # sza <- stats::na.omit(sza)

  # saa
  distrib_sun_angle$azimuth <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Azimuth$Values_List
  saa <- c()
  for (i in seq_len(length(distrib_sun_angle$zenith)))
    saa <- c(saa,as.numeric(strsplit(x = distrib_sun_angle$azimuth[i]$VALUES,
                                     split = ' ')[[1]]))
  # saa <- stats::na.omit(saa)

  # vza
  vza <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
      for (k in seq_len(length(values)))
        vza <- c(vza,as.numeric(strsplit(x = values[k]$VALUES,
                                         split = ' ')[[1]]))
    }
  }
  # vza <- stats::na.omit(vza)

  # vaa
  vaa <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
      for (k in seq_len(length(values)))
        vaa <- c(vaa,as.numeric(strsplit(x = values[k]$VALUES,
                                         split = ' ')[[1]]))
    }
  }
  # vaa <- stats::na.omit(vaa)
  return(list('saa' = saa, 'sza' = sza, 'vaa' = vaa, 'vza' = vza))
}
