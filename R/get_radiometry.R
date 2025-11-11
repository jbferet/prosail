#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param sensor_name character. name of the sensor
#' @param spectral_properties list. list of spectral properties including
#' wl = central wavelength and fwhm = corresponding fwhm
#' @param srf_path character. path for the file where srf should be
#' saved as CSV
#' @param save_srf boolean. Should srf be saved in file?
#'
#' @return srf list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @importFrom utils read.csv write.table
#' @export

get_spectral_response_function <- function(sensor_name = 'user_defined',
                           spectral_properties = NULL,
                           srf_path = './',
                           save_srf = TRUE){

  # == == == == == == == == == == == == == == == == == == == == == == == == =
  ### if the spectral response function of the sensor is already defined  ###
  # == == == == == == == == == == == == == == == == == == == == == == == == =
  if (toupper(sensor_name) %in% c('S2', 'SENTINEL2',
                                  'SENTINEL-2', 'SENTINEL_2')){
    srf <- prosail::Sentinel_2
    srf$sensor <- 'Sentinel_2'
  } else if (toupper(sensor_name) %in% c('S2A', 'SENTINEL2A',
                                         'SENTINEL-2A', 'SENTINEL_2A')){
    srf <- prosail::Sentinel_2A
    srf$sensor <- 'Sentinel_2A'
  } else if (toupper(sensor_name) %in% c('S2B', 'SENTINEL2B',
                                         'SENTINEL-2B', 'SENTINEL_2B')){
    srf <- prosail::Sentinel_2B
    srf$sensor <- 'Sentinel_2B'
  } else if (toupper(sensor_name) %in% c('S2C', 'SENTINEL2C',
                                         'SENTINEL-2C', 'SENTINEL_2C')){
    srf <- prosail::Sentinel_2C
    srf$sensor <- 'Sentinel_2C'
  } else if (toupper(sensor_name) =='VENUS'){
    srf <- prosail::Venus
    srf$sensor <- 'Venus'
  } else if (toupper(sensor_name) %in% c('L7', 'LANDSAT_7', 'LANDSAT-7')){
    srf <- prosail::Landsat_7
    srf$sensor <- 'Landsat_7'
  } else if (toupper(sensor_name) %in% c('L8', 'LANDSAT_8', 'LANDSAT-8')){
    srf <- prosail::Landsat_8
    srf$sensor <- 'Landsat_8'
  } else if (toupper(sensor_name) %in% c('L9', 'LANDSAT_9', 'LANDSAT-9')){
    srf <- prosail::Landsat_9
    srf$sensor <- 'Landsat_9'
  } else if (toupper(sensor_name) == 'MODIS'){
    srf <- prosail::MODIS
    srf$sensor <- 'MODIS'
  } else if (toupper(sensor_name) %in% c('SPOT_6_7','SPOT_6', 'SPOT6', 'SPOT-6',
                                         'SPOT_7', 'SPOT7', 'SPOT-7')){
    srf <- prosail::SPOT_6_7
    srf$sensor <- 'SPOT_6_7'
  } else if (toupper(sensor_name) %in% c('PLEIADES', 'PLEIADES_1A')){
    srf <- prosail::Pleiades_1A
    srf$sensor <- 'Pleiades_1A'
  } else if (toupper(sensor_name)=='PLEIADES_1B'){
    srf <- prosail::Pleiades_1B
    srf$sensor <- 'Pleiades_1B'
    # == == == == == == == == == == == == == == == == == == == == == == == == =
    ### if the spectral response function of the sensor is not defined      ###
    ###       but spectral characteristics (wl & fwhm) are provided         ###
    # == == == == == == == == == == == == == == == == == == == == == == == == =
    # if sensor is none of the above
  }  else {
    # check if spectral properties are provided in order to compute spectral
    # response based on gaussian assumption
    if (!is.null(spectral_properties)){
      if (!'wl' %in% names(spectral_properties) | ! 'fwhm' %in% names(spectral_properties)){
        message('Input variable "spectral_properties" expected to include fields')
        message('- wl (central wavelengths of bands in nm) ')
        message('- fwhl (fwhm of spectral bands in nm)')
        message('Please correct this input variable & provide proper info ')
        stop()
      } else if ('wl' %in% names(spectral_properties) &
                 'fwhm' %in% names(spectral_properties)){
        # check if proper spectral range
        if (max(spectral_properties$wl)<400 | min(spectral_properties$wl)>2500){
          message('Spectral bands in "spectral_properties" should be defined in nm')
          message('and should include bands between 400 and 2500')
          stop()
        }
        # produce spectral response per channel
        srf <- compute_srf(wvl = spectral_properties$wl,
                           fwhm = spectral_properties$fwhm,
                           sensor_name = sensor_name)
        # save srf as csv
        srf_save <- cbind(srf$original_bands,format(srf$spectral_response,
                                                   digits = 4,
                                                   scientific = FALSE))
        colnames(srf_save) <- c('SR_WL',srf$spectral_bands)
        if (save_srf==TRUE){
          path_srf <- file.path(srf_path,
                                paste0(sensor_name,'_spectral_response.csv'))
          message('Saving sensor spectral response function here:')
          print(path_srf)
          write.table(x = srf_save,file = path_srf,
                      quote = FALSE, sep = '\t',
                      col.names = TRUE, row.names = FALSE)
        }
      }
      srf$spectral_response <- t(srf$spectral_response)

      ## == == == == == == == == == == == == == == == == == == == == == == ==##
      ## if the spectral response function of the sensor is not defined      ##
      ##     and no spectral characteristics (wl & fwhm) are provided        ##
      ##   but a file corresponding to sensor characteristics is found       ##
      ## == == == == == == == == == == == == == == == == == == == == == == ==##
    } else {
      path_srf <- file.path(srf_path,
                            paste0(sensor_name,'_spectral_response.csv'))
      if (file.exists(path_srf)){
        # read file containing spectral response
        message('_____ reading spectral response corresponding to ______')
        print(sensor_name)
        srfraw <- read.csv(path_srf,header = TRUE,sep = '\t')
        original_bands <- srfraw[,1]
        spectral_bands <- colnames(srfraw)
        spectral_bands <- spectral_bands[-1]
        # check if conversion of spectral bands into numeric values
        spectral_bands_tmp <- as.numeric(gsub(pattern = 'X',
                                              replacement = '',
                                              x = spectral_bands))
        if (length(which(is.na(spectral_bands_tmp)))==0)
          spectral_bands <- spectral_bands_tmp
        spectral_response <- srfraw[,-1]
        spectral_response <- t(spectral_response)
        srf <- list('spectral_response' = spectral_response,
                    'spectral_bands' = spectral_bands,
                    'original_bands'= original_bands,
                    'central_wl' = spectral_bands,
                    'sensor' = sensor_name)

        ### == == == == == == == == == == == == == == == == == == == == == ==###
        ### if the spectral response function of the sensor is not defined   ###
        ###     and no spectral characteristics (wl & fwhm) are provided     ###
        ###   and no file corresponding to sensor characteristics is found   ###
        ### == == == == == == == == == == == == == == == == == == == == == ==###
      } else {
        message('___ Spectral response of the sensor expected here: ____')
        print(path_srf)
        message('_ Please check how to define sensor spectral response _')
        message('_____________ and create file accordingly _____________')
        srf <- NULL
      }
    }
  }
  return(srf)
}


#' @rdname prosail-deprecated
#' @export
GetRadiometry <- function(SensorName = 'Custom',
                          SpectralProps = NULL,
                          Path_SensorResponse = './',
                          SaveSRF = TRUE){
  .Deprecated("get_spectral_response_function")
  get_spectral_response_function(SensorName, SpectralProps, Path_SensorResponse, SaveSRF)
}
