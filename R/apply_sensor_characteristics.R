#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param wvl numeric. spectral sampling of the input spectral data
#' @param refl numeric. input spectral data (unique sample)
#' @param srf list. Spectral response & Spectral Bands of the sensor & Original
#' Bands for which srf is defined
#' @return refl_sensor numeric. Output spectral data, sensor resolution
#' @export

apply_sensor_characteristics <- function(wvl, refl, srf){

  refl <- data.frame(refl)
  nb_bands_origin <- length(wvl)
  if (dim(srf$spectral_response)[1]==nb_bands_origin)
    srf$spectral_response <- t(srf$spectral_response)
  nb_bands_sensor <- dim(srf$spectral_response)[1]
  array_input <- cbind(matrix(wvl,ncol = 1),refl)
  refl_sensor <- list()
  for (i in seq_len(nb_bands_sensor)){
    # for which spectral bands spectral response >0
    index_in <- which(srf$spectral_response[i,]>0)
    usfl_wvl <- srf$original_bands[index_in]
    # which spectral bands from refl correspond to sensor bands
    index_out <- which(is.element(wvl,usfl_wvl))
    # for the given band defined in index_in
    band_values <-srf$spectral_response[i,index_in]
    # simulated reflectance based on weighting
    if (length(index_in) == length(index_out)){
      pond <- matrix(band_values,
                     nrow = 1)%*%as.matrix(refl[index_out,])
      integrate_channel <- sum(band_values)
      refl_channel <- pond/integrate_channel
      refl_sensor[[i]] <- matrix(refl_channel,nrow =1)
    } else {
      message('sensor SRF should be compatible with prosail ')
      message(' (1 nm spectral sampling from 400 to 2500 m) ')
      message('         This is currently not the case      ')
      message('The follwing spectral band is out of spectral range')
      print(srf$spectral_bands[i])
      message('           Values will be set to 0           ')
      refl_sensor[[i]] <- matrix(0, ncol = ncol(refl), nrow =1)
    }
  }
  refl_sensor <-data.frame(do.call('rbind',refl_sensor))
  if (is.null(colnames(refl)))
    colnames(refl) <- paste0('sample#',
                                         seq_len(ncol(refl)))
  names(refl_sensor) <- colnames(refl)
  return(refl_sensor)
}


#' @rdname prosail-deprecated
#' @export
applySensorCharacteristics <- function(wvl, InRefl, SRF){
  .Deprecated("apply_sensor_characteristics")
  apply_sensor_characteristics(wvl, InRefl, SRF)
}
