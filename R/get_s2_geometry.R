#' This function returns geometry of acquisition for S2 image
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#'
#' @param MTD_TL_xml character. Path for metadata file MTD_TL.xml
#' @param verbose Boolean. Should messages be displayed?
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @importFrom XML xml xmlToList
#' @export

get_s2_geometry <- function(MTD_TL_xml, verbose = FALSE){

  # read XML file containing info about geometry of acquisition
  s2xml <- XML::xml(MTD_TL_xml)
  s2xml <- XML::xmlToList(s2xml)
  if (is.null(s2xml$Dataset_Identification$AUTHORITY)){
    geom_s2 <- get_s2_geometry_from_SAFE(s2xml)
  } else if (s2xml$Dataset_Identification$AUTHORITY=='THEIA'){
    if (verbose==TRUE){
      message('identification of S2 image produced by THEIA')
      message(s2xml$Dataset_Identification$IDENTIFIER)
    }
    geom_s2 <- get_s2_geometry_from_THEIA(s2xml)
  } else {
    geom_s2 <- get_s2_geometry_from_SAFE(s2xml)
  }
  return(list('SAA' = geom_s2$saa, 'SZA' = geom_s2$sza,
              'VAA' = geom_s2$vaa, 'VZA' = geom_s2$vza))
}
