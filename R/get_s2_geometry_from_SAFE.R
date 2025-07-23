#' This function returns geometry of acquisition for S2 image processed with
#' Sen2Cor
#' - sza = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - vza = list of viewer zenith angle
#' - vaa = list of viewer azimuth angle
#' @param s2xml list. produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (sza, saa, vza, vaa)
#' @importFrom stats na.omit
#' @export

get_s2_geometry_from_SAFE <- function(s2xml){

  distrib_sun_angle <- list()
  # sza
  distrib_sun_angle$zenith <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Zenith$Values_List
  sza <- list()
  for (i in seq_len(length(distrib_sun_angle$zenith))){
    sza[[i]] <- as.numeric(strsplit(x = distrib_sun_angle$zenith[i]$VALUES,
                                    split = ' ')[[1]])
  }
  sza <- do.call(what = rbind, args = sza)
  # sza <- stats::na.omit(sza)

  # saa
  distrib_sun_angle$azimuth <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Azimuth$Values_List
  saa <- list()
  for (i in seq_len(length(distrib_sun_angle$azimuth))){
    saa[[i]] <- as.numeric(strsplit(x = distrib_sun_angle$azimuth[i]$VALUES,
                                    split = ' ')[[1]])
  }
  saa <- do.call(what = rbind, args = saa)

  # vza
  vza <- list()
  ii <- 0
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    ii <- ii + 1
    vza_d <- band_detector[i]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
    for (j in seq_len(length(distrib_sun_angle$zenith))){
      vza_d[[j]] <- as.numeric(strsplit(x = vza_d[[j]],
                                      split = ' ')[[1]])
    }
    vza[[ii]] <- do.call(what = rbind, args = vza_d)
  }
  # get mean view zenith angle for all detectors
  vza <- apply(simplify2array(vza), 1:2, mean, na.rm = TRUE)

  # vaa
  vaa <- list()
  ii <- 0
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    ii <- ii + 1
    vaa_d <- band_detector[i]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
    for (j in seq_len(length(distrib_sun_angle$azimuth))){
      vaa_d[[j]] <- as.numeric(strsplit(x = vaa_d[[j]],
                                        split = ' ')[[1]])
    }
    vaa[[ii]] <- do.call(what = rbind, args = vaa_d)
  }
  # get mean view azimuth angle for all detectors
  vaa <- apply(simplify2array(vaa), 1:2, mean, na.rm = TRUE)
  return(list('saa' = saa, 'sza' = sza, 'vaa' = vaa, 'vza' = vza))
}


#' @rdname prosail-deprecated
#' @export
get_S2geometry_from_SAFE <- function(s2xml){
  .Deprecated("get_s2_geometry_from_SAFE")
  get_s2_geometry_from_SAFE(s2xml)
}
