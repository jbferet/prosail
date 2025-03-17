#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param SensorName character. name of the sensor
#' @param SpectralProps list. list of spectral properties including
#' wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character. path for the file where SRF should be
#' saved as CSV
#' @param SaveSRF boolean. Should SRF be saved in file?
#'
#' @return SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @importFrom utils read.csv write.table
#' @export

get_radiometry <- function(SensorName = 'Custom',
                          SpectralProps = NULL,
                          Path_SensorResponse = './',
                          SaveSRF = TRUE){

  # == == == == == == == == == == == == == == == == == == == == == == == == =
  ### if the spectral response function of the sensor is already defined  ###
  # == == == == == == == == == == == == == == == == == == == == == == == == =
  if (toupper(SensorName) %in% c('S2', 'SENTINEL2',
                                 'SENTINEL-2', 'SENTINEL_2')){
    SRF <- prosail::Sentinel_2
    SRF$Sensor <- 'Sentinel_2'
  } else if (toupper(SensorName) %in% c('S2A', 'SENTINEL2A',
                                        'SENTINEL-2A', 'SENTINEL_2A')){
    SRF <- prosail::Sentinel_2A
    SRF$Sensor <- 'Sentinel_2A'
  } else if (toupper(SensorName) %in% c('S2B', 'SENTINEL2B',
                                        'SENTINEL-2B', 'SENTINEL_2B')){
    SRF <- prosail::Sentinel_2B
    SRF$Sensor <- 'Sentinel_2B'
  } else if (toupper(SensorName) =='VENUS'){
    SRF <- prosail::Venus
    SRF$Sensor <- 'Venus'
  } else if (toupper(SensorName) %in% c('L7', 'LANDSAT_7', 'LANDSAT-7')){
    SRF <- prosail::Landsat_7
    SRF$Sensor <- 'Landsat_7'
  } else if (toupper(SensorName) %in% c('L8', 'LANDSAT_8', 'LANDSAT-8')){
    SRF <- prosail::Landsat_8
    SRF$Sensor <- 'Landsat_8'
  } else if (toupper(SensorName) %in% c('L9', 'LANDSAT_9', 'LANDSAT-9')){
    SRF <- prosail::Landsat_9
    SRF$Sensor <- 'Landsat_9'
  } else if (toupper(SensorName) == 'MODIS'){
    SRF <- prosail::MODIS
    SRF$Sensor <- 'MODIS'
  } else if (toupper(SensorName) %in% c('SPOT_6_7', 'SPOT_6', 'SPOT6', 'SPOT-6',
                                        'SPOT_7', 'SPOT7', 'SPOT-7')){
    SRF <- prosail::SPOT_6_7
    SRF$Sensor <- 'SPOT_6_7'
  } else if (toupper(SensorName) %in% c('PLEIADES', 'PLEIADES_1A')){
    SRF <- prosail::Pleiades_1A
    SRF$Sensor <- 'Pleiades_1A'
  } else if (toupper(SensorName)=='PLEIADES_1B'){
    SRF <- prosail::Pleiades_1B
    SRF$Sensor <- 'Pleiades_1B'
    # == == == == == == == == == == == == == == == == == == == == == == == == =
    ### if the spectral response function of the sensor is not defined      ###
    ###       but spectral characteristics (wl & FWHM) are provided         ###
    # == == == == == == == == == == == == == == == == == == == == == == == == =
    # if sensor is none of the above
  }  else {
    # check if spectral properties are provided in order to compute spectral
    # response based on gaussian assumption
    if (!is.null(SpectralProps)){
      if (!'wl' %in% names(SpectralProps) | ! 'fwhm' %in% names(SpectralProps)){
        message('Input variable "SpectralProps" expected to include fields')
        message('- wl (central wavelengths of bands in nm) ')
        message('- fwhl (FWHM of spectral bands in nm)')
        message('Please correct this input variable & provide proper info ')
        stop()
      } else if ('wl' %in% names(SpectralProps) &
                 'fwhm' %in% names(SpectralProps)){
        # check if proper spectral range
        if (max(SpectralProps$wl)<400 | min(SpectralProps$wl)>2500){
          message('Spectral bands in "SpectralProps" should be defined in nm')
          message('and should include bands between 400 and 2500')
          stop()
        }
        # produce spectral response per channel
        SRF <- compute_SRF(wvl = SpectralProps$wl,
                           FWHM = SpectralProps$fwhm,
                           SensorName = SensorName)
        # save SRF as csv
        SRFsave <- cbind(SRF$Original_Bands,format(SRF$Spectral_Response,
                                                   digits = 4,
                                                   scientific = FALSE))
        colnames(SRFsave) <- c('SR_WL',SRF$Spectral_Bands)
        if (SaveSRF==TRUE){
          Path_SRF <- file.path(Path_SensorResponse,
                                paste0(SensorName,'_Spectral_Response.csv'))
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
      ##     and no spectral characteristics (wl & FWHM) are provided        ##
      ##   but a file corresponding to sensor characteristics is found       ##
      ## == == == == == == == == == == == == == == == == == == == == == == ==##
    } else {
      Path_SRF <- file.path(Path_SensorResponse,
                            paste0(SensorName,'_Spectral_Response.csv'))
      if (file.exists(Path_SRF)){
        # read file containing spectral response
        message('_____ reading spectral response corresponding to ______')
        print(SensorName)
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
                    'Sensor' = SensorName)

        ### == == == == == == == == == == == == == == == == == == == == == ==###
        ### if the spectral response function of the sensor is not defined   ###
        ###     and no spectral characteristics (wl & FWHM) are provided     ###
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
