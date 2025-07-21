#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param wvl numeric. spectral sampling of the input spectral data
#' @param input_refl_table numeric. input spectral data (unique sample)
#' @param srf list. Spectral response & Spectral Bands of the sensor & Original
#' Bands for which srf is defined
#' @return refl_sensor numeric. Output spectral data, sensor resolution
#' @export

apply_sensor_characteristics <- function(wvl, input_refl_table, srf){

  input_refl_table <- data.frame(input_refl_table)
  nb_bands_origin <- length(wvl)
  if (dim(srf$Spectral_Response)[1]==nb_bands_origin)
    srf$Spectral_Response <- t(srf$Spectral_Response)
  nb_bands_sensor <- dim(srf$Spectral_Response)[1]
  array_input <- cbind(matrix(wvl,ncol = 1),input_refl_table)
  refl_sensor <- list()
  for (i in seq_len(nb_bands_sensor)){
    # for which spectral bands spectral response >0
    index_in <- which(srf$Spectral_Response[i,]>0)
    usfl_wvl <- srf$Original_Bands[index_in]
    # which spectral bands from input_refl_table correspond to sensor bands
    index_out <- which(is.element(wvl,usfl_wvl))
    # for the given band defined in index_in
    band_values <-srf$Spectral_Response[i,index_in]
    # simulated reflectance based on weighting
    if (length(index_in) == length(index_out)){
      pond <- matrix(band_values,
                     nrow = 1)%*%as.matrix(input_refl_table[index_out,])
      integrate_channel <- sum(band_values)
      refl_channel <- pond/integrate_channel
      refl_sensor[[i]] <- matrix(refl_channel,nrow =1)
    } else {
      message('sensor SRF should be compatible with prosail ')
      message(' (1 nm spectral sampling from 400 to 2500 m) ')
      message('         This is currently not the case      ')
      message('The follwing spectral band is out of spectral range')
      print(srf$Spectral_Bands[i])
      message('           Values will be set to 0           ')
      refl_sensor[[i]] <- matrix(0, ncol = ncol(input_refl_table), nrow =1)
    }
  }
  refl_sensor <-data.frame(do.call('rbind',refl_sensor))
  if (is.null(colnames(input_refl_table)))
    colnames(input_refl_table) <- paste0('sample#',
                                         seq_len(ncol(input_refl_table)))
  names(refl_sensor) <- colnames(input_refl_table)
  return(refl_sensor)
}
