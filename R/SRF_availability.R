#' This function returns the list of al sensors for which SRF is available
#'
#' @return list of sensors
#' @export

srf_availability <- function(){
  return(c('SENTINEL-2', 'SENTINEL-2A', 'SENTINEL-2B', 'SENTINEL-2C', 'VENUS',
           'LANDSAT-7', 'LANDSAT-8', 'LANDSAT-9', 'MODIS',
           'SPOT-6', 'SPOT-7', 'PLEIADES', 'PLEIADES-1A', 'PLEIADES-1B'))
}
