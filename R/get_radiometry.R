#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param sensor_name character. name of the sensor
#' @param spectral_properties list. list of spectral properties including
#' wl = central wavelength and fwhm = corresponding fwhm
#' @param srf_path character. path for the file where SRF should be
#' saved as CSV
#' @param save_srf boolean. Should SRF be saved in file?
#'
#' @return SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @importFrom utils read.csv write.table
#' @export

get_radiometry <- function(sensor_name = 'Custom',
                           spectral_properties = NULL,
                           srf_path = './',
                           save_srf = TRUE){

  # == == == == == == == == == == == == == == == == == == == == == == == == =
  ### if the spectral response function of the sensor is already defined  ###
  # == == == == == == == == == == == == == == == == == == == == == == == == =
  if (toupper(sensor_name) %in% c('S2', 'SENTINEL2',
                                  'SENTINEL-2', 'SENTINEL_2')){
    SRF <- prosail::Sentinel_2
    SRF$Sensor <- 'Sentinel_2'
  } else if (toupper(sensor_name) %in% c('S2A', 'SENTINEL2A',
                                         'SENTINEL-2A', 'SENTINEL_2A')){
    SRF <- prosail::Sentinel_2A
    SRF$Sensor <- 'Sentinel_2A'
  } else if (toupper(sensor_name) %in% c('S2B', 'SENTINEL2B',
                                         'SENTINEL-2B', 'SENTINEL_2B')){
    SRF <- prosail::Sentinel_2B
    SRF$Sensor <- 'Sentinel_2B'
  } else if (toupper(sensor_name) %in% c('S2C', 'SENTINEL2C',
                                         'SENTINEL-2C', 'SENTINEL_2C')){
    SRF <- prosail::Sentinel_2C
    SRF$Sensor <- 'Sentinel_2C'
  } else if (toupper(sensor_name) =='VENUS'){
    SRF <- prosail::Venus
    SRF$Sensor <- 'Venus'
  } else if (toupper(sensor_name) %in% c('L7', 'LANDSAT_7', 'LANDSAT-7')){
    SRF <- prosail::Landsat_7
    SRF$Sensor <- 'Landsat_7'
  } else if (toupper(sensor_name) %in% c('L8', 'LANDSAT_8', 'LANDSAT-8')){
    SRF <- prosail::Landsat_8
    SRF$Sensor <- 'Landsat_8'
  } else if (toupper(sensor_name) %in% c('L9', 'LANDSAT_9', 'LANDSAT-9')){
    SRF <- prosail::Landsat_9
    SRF$Sensor <- 'Landsat_9'
  } else if (toupper(sensor_name) == 'MODIS'){
    SRF <- prosail::MODIS
    SRF$Sensor <- 'MODIS'
  } else if (toupper(sensor_name) %in% c('SPOT_6_7','SPOT_6', 'SPOT6', 'SPOT-6',
                                         'SPOT_7', 'SPOT7', 'SPOT-7')){
    SRF <- prosail::SPOT_6_7
    SRF$Sensor <- 'SPOT_6_7'
  } else if (toupper(sensor_name) %in% c('PLEIADES', 'PLEIADES_1A')){
    SRF <- prosail::Pleiades_1A
    SRF$Sensor <- 'Pleiades_1A'
  } else if (toupper(sensor_name)=='PLEIADES_1B'){
    SRF <- prosail::Pleiades_1B
    SRF$Sensor <- 'Pleiades_1B'
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
        SRF <- compute_SRF(wvl = spectral_properties$wl,
                           fwhm = spectral_properties$fwhm,
                           sensor_name = sensor_name)
        # save SRF as csv
        SRFsave <- cbind(SRF$Original_Bands,format(SRF$Spectral_Response,
                                                   digits = 4,
                                                   scientific = FALSE))
        colnames(SRFsave) <- c('SR_WL',SRF$Spectral_Bands)
        if (save_srf==TRUE){
          Path_SRF <- file.path(srf_path,
                                paste0(sensor_name,'_Spectral_Response.csv'))
          message('Saving sensor spectral response function here:')
          print(Path_SRF)
          write.table(x = SRFsave,file = Path_SRF,
                      quote = FALSE, sep = '\t',
                      col.names = TRUE, row.names = FALSE)
        }
      }
      SRF$Spectral_Response <- t(SRF$Spectral_Response)

      ## == == == == == == == == == == == == == == == == == == == == == == ==##
      ## if the spectral response function of the sensor is not defined      ##
      ##     and no spectral characteristics (wl & fwhm) are provided        ##
      ##   but a file corresponding to sensor characteristics is found       ##
      ## == == == == == == == == == == == == == == == == == == == == == == ==##
    } else {
      Path_SRF <- file.path(srf_path,
                            paste0(sensor_name,'_Spectral_Response.csv'))
      if (file.exists(Path_SRF)){
        # read file containing spectral response
        message('_____ reading spectral response corresponding to ______')
        print(sensor_name)
        SRFraw <- read.csv(Path_SRF,header = TRUE,sep = '\t')
        Original_Bands <- SRFraw[,1]
        Spectral_Bands <- colnames(SRFraw)
        Spectral_Bands <- Spectral_Bands[-1]
        # check if conversion of spectral bands into numeric values
        Spectral_Bands_tmp <- as.numeric(gsub(pattern = 'X',
                                              replacement = '',
                                              x = Spectral_Bands))
        if (length(which(is.na(Spectral_Bands_tmp)))==0)
          Spectral_Bands <- Spectral_Bands_tmp
        Spectral_Response <- SRFraw[,-1]
        Spectral_Response <- t(Spectral_Response)
        SRF <- list('Spectral_Response' = Spectral_Response,
                    'Spectral_Bands' = Spectral_Bands,
                    'Original_Bands'= Original_Bands,
                    'Central_WL' = Spectral_Bands,
                    'Sensor' = sensor_name)

        ### == == == == == == == == == == == == == == == == == == == == == ==###
        ### if the spectral response function of the sensor is not defined   ###
        ###     and no spectral characteristics (wl & fwhm) are provided     ###
        ###   and no file corresponding to sensor characteristics is found   ###
        ### == == == == == == == == == == == == == == == == == == == == == ==###
      } else {
        message('___ Spectral response of the sensor expected here: ____')
        print(Path_SRF)
        message('_ Please check how to define sensor spectral response _')
        message('_____________ and create file accordingly _____________')
        SRF <- NULL
      }
    }
  }
  return(SRF)
}
