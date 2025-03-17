#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param wvl numeric. spectral sampling of the input spectral data
#' @param InRefl numeric. input spectral data (unique sample)
#' @param SRF list. Spectral response & Spectral Bands of the sensor & Original
#' Bands for which SRF is defined
#' @return OutRefl numeric. Output spectral data, sensor resolution
#' @export

apply_sensor_characteristics <- function(wvl, InRefl, SRF){

  InRefl <- data.frame(InRefl)
  nbBands_Origin <- length(wvl)
  if (dim(SRF$Spectral_Response)[1]==nbBands_Origin){
    SRF$Spectral_Response <- t(SRF$Spectral_Response)
  }
  nbBands_Sensor <- dim(SRF$Spectral_Response)[1]
  array_input <- cbind(matrix(wvl,ncol = 1),InRefl)
  OutRefl <- list()
  for (i in seq_len(nbBands_Sensor)){
    # for which spectral bands spectral response >0
    indexIN <- which(SRF$Spectral_Response[i,]>0)
    usfl_wvl <- SRF$Original_Bands[indexIN]
    # which spectral bands from InRefl correspond to sensor spectral bands
    indexOUT <- which(is.element(wvl,usfl_wvl))
    # for the given band defined in indexIN
    band_values <-SRF$Spectral_Response[i,indexIN]
    # simulated reflectance based on weighting
    if (length(indexIN) == length(indexOUT)){
      pond <- matrix(band_values,nrow = 1)%*%as.matrix(InRefl[indexOUT,])
      IntegrateChannel <- sum(band_values)
      Refl_Channel <- pond/IntegrateChannel
      OutRefl[[i]] <- matrix(Refl_Channel,nrow =1)
    } else {
      message('Please make sure that the spectral resolution for sensor response function ')
      message('__ is compatible ewith model (1 nm spectral sampling from 400 to 2500 m) __')
      message('______________________ This is currently not the case _____________________')
      message('_____ The follwing spectral band is not within accepted spectral range ____')
      print(SRF$Spectral_Bands[i])
      message('_________________________ Values will be set to 0  ________________________')
      OutRefl[[i]] <- matrix(0, ncol = ncol(InRefl), nrow =1)
    }
  }
  OutRefl <-data.frame(do.call('rbind',OutRefl))
  if (is.null(colnames(InRefl)))
    colnames(InRefl) <- paste0('sample#',seq_len(ncol(InRefl)))
  names(OutRefl) <- colnames(InRefl)
  return(OutRefl)
}
