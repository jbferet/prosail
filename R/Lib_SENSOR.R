# ============================================================================= =
# prosail
# Lib_SENSOR.R
# ============================================================================= =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================= =
# This Library includes functions dedicated to conversion from high resolution
# reflectance to sensor reflectance using either known sensor response,
# or gaussian filters
# ============================================================================= =

#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param SensorName character. name of the sensor
#' @param SpectralProps list. list of spectral properties including
#' wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character. path for the file where SRF should be
#' saved as CSV
#' @param SaveSRF boolean. Should SRF be saved in file?
#' @return SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
# @import utils
#' @importFrom utils read.csv write.table
#' @export

GetRadiometry <- function(SensorName = 'Custom',
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
  } else if (toupper(SensorName) %in% c('S2C', 'SENTINEL2C',
                                        'SENTINEL-2C', 'SENTINEL_2C')){
    SRF <- prosail::Sentinel_2C
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
        SRF <- Compute_SRF(wvl = SpectralProps$wl,
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

      ### == == == == == == == == == == == == == == == == == == == == == == ==###
      ### if the spectral response function of the sensor is not defined      ###
      ###     and no spectral characteristics (wl & FWHM) are provided        ###
      ###   but a file corresponding to sensor characteristics is found       ###
      ### == == == == == == == == == == == == == == == == == == == == == == ==###
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
        if (length(which(is.na(Spectral_Bands_tmp)))==0) Spectral_Bands <- Spectral_Bands_tmp
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

#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param wvl numeric. spectral sampling of the input spectral data
#' @param InRefl numeric. input spectral data (unique sample)
#' @param SRF list. Spectral response & Spectral Bands of the sensor & Original
#' Bands for which SRF is defined
#' @return OutRefl numeric. Output spectral data, sensor resolution
#' @export

applySensorCharacteristics <- function(wvl, InRefl, SRF){

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
  if (is.null(colnames(InRefl))) colnames(InRefl) <- paste0('sample#',seq_len(ncol(InRefl)))
  names(OutRefl) <- colnames(InRefl)
  return(OutRefl)
}

#' This function converts high resolution spectral info into broader
#' spectral characteristics
#'
#' @param SpecPROSPECT list. optical constants required for PROSPECT
#' @param SpecSOIL list. dry soil and wet soil, or a unique soil sample if
#'  psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
#' @param SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#'
#' @return SpecSensor list. list of input specs with sensor resolution
#' @export

PrepareSensorSimulation <- function(SpecPROSPECT,SpecSOIL,SpecATM,SRF){

  # adjust optical constants
  wvl <- SpecPROSPECT$lambda
  # leaf properties
  SpecPROSPECT_Sensor <- applySensorCharacteristics(wvl,SpecPROSPECT,SRF)
  # atmospheric properties
  SpecATM_Sensor <- applySensorCharacteristics(wvl,SpecATM,SRF)
  # soil properties
  SpecSOIL_Sensor <- applySensorCharacteristics(wvl,SpecSOIL,SRF)
  SpecSensor <- list('SpecPROSPECT_Sensor' = SpecPROSPECT_Sensor,
                     'SpecATM_Sensor' = SpecATM_Sensor,
                     'SpecSOIL_Sensor' = SpecSOIL_Sensor,
                     'BandNames' = SRF$Spectral_Bands)
  return(SpecSensor)
}

#' Computes spectral response function based on wavelength and FWHM characteristics
#' @param wvl numeric. spectral sampling of the sensor
#' @param FWHM numeric. Full Width Half Maximum for each spectral band
#' @param SensorName character. sensor name
#'
#' @return SRF list. spectral response, sensor Spectral Bands & Original Bands
#' for which SRF is defined
#' @importFrom stats dnorm
#' @export

Compute_SRF <- function(wvl,FWHM, SensorName = 'Custom'){

  # define full spectral domain in optical domain
  lambda <- seq(400,2500,by = 1)
  Spectral_Response <- matrix(0,ncol = length(wvl),nrow = length(lambda))
  for (i in seq_len(length(wvl))){
    y <- dnorm(lambda,wvl[i],FWHM[i]/2.355)
    Spectral_Response[,i] <- y/max(y)
  }
  Spectral_Response[which(Spectral_Response < 0.001)] <- 0
  SRF <- list('Spectral_Response' = Spectral_Response,
              'Spectral_Bands' = wvl,
              'Original_Bands' = lambda,
              'Central_WL' = wvl,
              'Sensor' = SensorName)
  return(SRF)
}

#' Computes spectral sensor characteristics for optical constants required for PROSPECT
#' soil parameters as well as direct and diffuse radiation for clear conditions
#' @param SensorName character. name of the sensor. Either already defined, or
#' to be defined based on other inputs
#' @param SpectralProps list. Sensor spectral characteristics including
#' wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character. Path where to get or save SRF
#'
#' @return SRF list. Spectral response function, corresponding spectral bands,
#' and Original Bands
#' @export

get_spec_sensor <- function(SensorName = 'MyCustomSensor',
                            SpectralProps = NULL,
                            Path_SensorResponse = NULL){
  SRF <- GetRadiometry(SensorName,SpectralProps = SpectralProps,
                       Path_SensorResponse = Path_SensorResponse)
  if (is.null(SRF)) stop()
  # adjust optical constants from 1nm sampling into spectral S2 spectral sampling
  SpecSensor <- PrepareSensorSimulation(prosail::SpecPROSPECT_FullRange,
                                        prosail::SpecSOIL,
                                        prosail::SpecATM,SRF)
  return(SpecSensor)
}

#' This function returns geometry of acquisition for S2 image
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#' @param MTD_TL_xml character. Path for metadata file MTD_TL.xml
#' @param verbose Boolean. Should messages be displayed?
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @importFrom XML xml xmlToList
#' @export

get_S2geometry <- function(MTD_TL_xml,verbose=FALSE){

  # read XML file containing info about geometry of acquisition
  s2xml <- XML::xml(MTD_TL_xml)
  s2xml <- XML::xmlToList(s2xml)
  if (is.null(s2xml$Dataset_Identification$AUTHORITY)){
    GeomS2 <- get_S2geometry_from_SAFE(s2xml)
  } else if (s2xml$Dataset_Identification$AUTHORITY=='THEIA'){
    if (verbose==TRUE){
      message('identification of S2 image produced by THEIA')
      message(s2xml$Dataset_Identification$IDENTIFIER)
    }
    GeomS2 <- get_S2geometry_from_THEIA(s2xml)
  } else {
    GeomS2 <- get_S2geometry_from_SAFE(s2xml)
  }
  return(list('SAA' = GeomS2$SAA, 'SZA' = GeomS2$SZA,
              'VAA' = GeomS2$VAA, 'VZA' = GeomS2$VZA))
}

#' This function returns geometry of acquisition for S2 image processed with MAJA
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#' @param s2xml list. list produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @export

get_S2geometry_from_THEIA <- function(s2xml){

  Distrib_SunAngle <- list()
  # SZA
  Distrib_SunAngle$Zenith <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Zenith$Values_List
  SZA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SZA <- c(SZA, as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,split = ' ')[[1]]))
  }
  SZA <- SZA[which(!is.na(SZA))]

  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Informations$Angles_Grids_List$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SAA <- c(SAA,as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,split = ' ')[[1]]))
  }
  SAA <- SAA[which(!is.na(SAA))]

  # VZA
  VZA <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
      for (k in seq_len(length(values))) {
        VZA <- c(VZA,as.numeric(strsplit(x = values[k]$VALUES,split = ' ')[[1]]))
      }
    }
  }
  VZA <- VZA[which(!is.na(VZA))]

  # VAA
  VAA <- c()
  band <- s2xml$Geometric_Informations$Angles_Grids_List$Viewing_Incidence_Angles_Grids
  for (i in seq_len(length(band))) {
    detector <- band[i]$Band_Viewing_Incidence_Angles_Grids_List
    for (j in seq_len((length(detector)-1))) {
      values  <- detector[j]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
      for (k in seq_len(length(values))) {
        VAA <- c(VAA,as.numeric(strsplit(x = values[k]$VALUES,split = ' ')[[1]]))
      }
    }
  }
  VAA <- VAA[which(!is.na(VAA))]
  return(list('SAA' = SAA, 'SZA' = SZA, 'VAA' = VAA, 'VZA' = VZA))
}

#' This function returns geometry of acquisition for S2 image processed with Sen2Cor
#' - SZA = list of sun zenith angle
#' - SAA = list of sun azimuth angle
#' - VZA = list of viewer zenith angle
#' - VAA = list of viewer azimuth angle
#' @param s2xml list. list produced from reading XML metadata file with package XML
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @export

get_S2geometry_from_SAFE <- function(s2xml){

  Distrib_SunAngle <- list()
  # SZA
  Distrib_SunAngle$Zenith <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Zenith$Values_List
  SZA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SZA <- c(SZA,as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,
                                     split = ' ')[[1]]))
  }
  SZA <- SZA[which(!is.na(SZA))]

  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- c()
  for (i in seq_len(length(Distrib_SunAngle$Zenith))){
    SAA <- c(SAA,as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,
                                     split = ' ')[[1]]))
  }
  SAA <- SAA[which(!is.na(SAA))]

  # VZA
  VZA <- c()
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    values <- band_detector[i]$Viewing_Incidence_Angles_Grids$Zenith$Values_List
    for (j in seq_len(length(values))){
      VZA <- c(VZA,as.numeric(strsplit(x = values[j]$VALUES, split = ' ')[[1]]))
    }
  }
  VZA <- VZA[which(!is.na(VZA))]

  # VAA
  VAA <- c()
  band_detector <- s2xml$Geometric_Info$Tile_Angles
  for (i in 3:(length(band_detector)-2)) {
    values <- band_detector[i]$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
    for (j in seq_len(length(values))){
      VAA <- c(VAA,as.numeric(strsplit(x = values[j]$VALUES,split = ' ')[[1]]))
    }
  }
  VAA <- VAA[which(!is.na(VAA))]
  return(list('SAA' = SAA, 'SZA' = SZA, 'VAA' = VAA, 'VZA' = VZA))
}

#' This function returns the list of al sensors for which SRF is available
#'
#' @return list of sensors
#' @export

SRF_Available <- function(){
  return(c('SENTINEL-2', 'SENTINEL-2A', 'SENTINEL-2B', 'VENUS',
           'LANDSAT-7', 'LANDSAT-8', 'LANDSAT-9', 'MODIS',
           'SPOT-6', 'SPOT-7', 'PLEIADES', 'PLEIADES-1A', 'PLEIADES-1B'))
}
