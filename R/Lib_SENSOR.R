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
#' @param SpectralProps list. list of spectral properties including wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character.
#' @return SRF list. Spectral response function, corresponding spectral bands, and Original Bands
# @import utils
#' @importFrom utils read.csv write.table
#' @export
GetRadiometry <- function(SensorName = 'MyCustomSensor',
                          SpectralProps = NULL,
                          Path_SensorResponse = NULL){

  # == == == == == == == == == == == == == == == == == == == == == == == == =
  ### if the spectral response function of the sensor is already defined  ###
  # == == == == == == == == == == == == == == == == == == == == == == == == =
  if (SensorName=='Sentinel_2'){                                  # if sensor is SENTINEL-2
    Spectral_Response <- prosail::Sentinel_2$Spectral_Response
    Spectral_Bands <- prosail::Sentinel_2$Spectral_Bands
    OriginalBands <- prosail::Sentinel_2$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='Sentinel_2A'){                          # if sensor is SENTINEL-2A
    Spectral_Response <- prosail::Sentinel_2A$Spectral_Response
    Spectral_Bands <- prosail::Sentinel_2A$Spectral_Bands
    OriginalBands <- prosail::Sentinel_2A$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='Sentinel_2B'){                          # if sensor is SENTINEL-2B
    Spectral_Response <- prosail::Sentinel_2B$Spectral_Response
    Spectral_Bands <- prosail::Sentinel_2B$Spectral_Bands
    OriginalBands <- prosail::Sentinel_2B$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='Venus'){                                # if sensor is Venus
    Spectral_Response <- prosail::Venus$Spectral_Response
    Spectral_Bands <- prosail::Venus$Spectral_Bands
    OriginalBands <- prosail::Venus$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='Landsat_7'){                            # if sensor is Landsat-7
    Spectral_Response <- prosail::Landsat_7$Spectral_Response
    Spectral_Bands <- prosail::Landsat_7$Spectral_Bands
    OriginalBands <- prosail::Landsat_7$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='Landsat_8'){                            # if sensor is Landsat-8
    Spectral_Response <- prosail::Landsat_8$Spectral_Response
    Spectral_Bands <- prosail::Landsat_8$Spectral_Bands
    OriginalBands <- prosail::Landsat_8$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='MODIS'){                                # if sensor is MODIS
    Spectral_Response <- prosail::MODIS$Spectral_Response
    Spectral_Bands <- prosail::MODIS$Spectral_Bands
    OriginalBands <- prosail::MODIS$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
  } else if (SensorName=='SPOT_6_7'){                             # if sensor is SPOT_6_7
    Spectral_Response <- prosail::SPOT_6_7$Spectral_Response
    Spectral_Bands <- prosail::SPOT_6_7$Spectral_Bands
    OriginalBands <- prosail::SPOT_6_7$OriginalBands
    SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)

    # == == == == == == == == == == == == == == == == == == == == == == == == =
    ### if the spectral response function of the sensor is not defined      ###
    ###       but spectral characteristics (wl & FWHM) are provided         ###
    # == == == == == == == == == == == == == == == == == == == == == == == == =
  }  else {                                                       # if sensor is none of the above
    # check if spectral properties are provided in order to compute spectral response based on gaussian assumption
    if (!is.null(SpectralProps)){
      if (is.na(match('wl',names(SpectralProps))) | is.na(match('fwhm',names(SpectralProps))>0 )){
        message('Input variable "SpectralProps" expected to include fields')
        message('- wl (central wavelengths of bands in nm) ')
        message('- fwhl (FWHM of spectral bands in nm)')
        message('Please correct this input variable & provide proper info ')
        stop()
      } else if (match('wl',names(SpectralProps))>0 & match('fwhm',names(SpectralProps))>0 ){
        # check if proper spectral range
        if (max(SpectralProps$wl)<400 | min(SpectralProps$wl)>2500){
          message('Spectral bands in "SpectralProps" should be defined in nm')
          message('and should include bands between 400 and 2500')
          stop()
        }
        # produce spectral response per channel
        SRF <- Compute_SRF(wvl = SpectralProps$wl,FWHM = SpectralProps$fwhm)
        # save SRF as csv
        SRFsave <- cbind(SRF$OriginalBands,format(SRF$Spectral_Response, digits = 4, scientific = FALSE))
        colnames(SRFsave) <- c('SR_WL',SRF$Spectral_Bands)
        if (!is.null(Path_SensorResponse)){                           # if a path is provided to find the sensor
          Path_SRF <- file.path(Path_SensorResponse,paste(SensorName,'_Spectral_Response.csv',sep=''))
          message('Saving spectral response function of sensor ')
          message('in following directory defined by "Path_SensorResponse" : ')
        } else {                                                      # if no path is provided to find the sensor, assuming it is in the WD
          Path_SRF <- file.path(paste(SensorName,'_Spectral_Response.csv',sep=''))
          message('Saving spectral response function of sensor ')
          message('in current working directory as file : ')
        }
        print(Path_SRF)
        write.table(x = SRFsave,file = Path_SRF,
                    quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
      }
      SRF$Spectral_Response <- t(SRF$Spectral_Response)

      ### == == == == == == == == == == == == == == == == == == == == == == ==###
      ### if the spectral response function of the sensor is not defined      ###
      ###     and no spectral characteristics (wl & FWHM) are provided        ###
      ###   but a file corresponding to sensor characteristics is found       ###
      ### == == == == == == == == == == == == == == == == == == == == == == ==###
    } else {
      if (!is.null(Path_SensorResponse)){                           # if a path is provided to find the sensor
        Path_SRF <- file.path(Path_SensorResponse,paste(SensorName,'_Spectral_Response.csv',sep=''))
      } else {                                                      # if no path is provided to find the sensor, assuming it is in the WD
        Path_SRF <- file.path(paste(SensorName,'_Spectral_Response.csv',sep=''))
      }
      if (file.exists(Path_SRF)){
        # read file containing spectral response
        message('_____ reading spectral response corresponding to ______')
        print(SensorName)
        SRFraw <- read.csv(Path_SRF,header = TRUE,sep = '\t')
        OriginalBands <- SRFraw[,1]
        Spectral_Bands <- colnames(SRFraw)
        Spectral_Bands <- Spectral_Bands[-1]
        # check if conversion of spctral bands into numeric values
        Spectral_Bands_tmp = as.numeric(gsub(pattern = 'X',replacement = '',x = Spectral_Bands))
        if (length(which(is.na(Spectral_Bands_tmp)))==0){
          Spectral_Bands = Spectral_Bands_tmp
        }
        Spectral_Response <-SRFraw[,-1]
        Spectral_Response <- t(Spectral_Response)
        SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)

        ### == == == == == == == == == == == == == == == == == == == == == == ==###
        ### if the spectral response function of the sensor is not defined      ###
        ###     and no spectral characteristics (wl & FWHM) are provided        ###
        ###   and no file corresponding to sensor characteristics is found      ###
        ### == == == == == == == == == == == == == == == == == == == == == == ==###
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
#' @param SRF list. Information about spectral response Spectral Bands of the sensor and Original Bands for which SRF is defined
#' @return OutRefl numeric. Output spectral data, sensor resolution
#' @export
applySensorCharacteristics <- function(wvl,InRefl,SRF){

  nbBands_Origin <- length(wvl)
  nbBands_Sensor <- dim(SRF$Spectral_Response)[1]
  # if (!length(InRefl)==nbBands_Origin){
  #   message('_ Dimensions of input data and spectral sampling do not match _')
  #   message('___________________ Please adjust input data___________________')
  #   message('_____ Input spectral data should be an individual sample ______')
  #   stop()
  # }
  array_input <- cbind(matrix(wvl,ncol = 1),InRefl)
  OutRefl <- list()
  for (i in 1:nbBands_Sensor){
    # identify on which spectral bands spectral response >0
    indexIN <- which(SRF$Spectral_Response[i,]>0)
    usfl_wvl <- SRF$OriginalBands[indexIN]
    # identify which spectral bands from InRefl correspond to the spectral bands from sensor
    indexOUT <- which(is.element(wvl,usfl_wvl))
    # for the given band defined in indexIN
    band_values <-SRF$Spectral_Response[i,indexIN]
    # simulated reflectance based on weighting
    if (length(indexIN)==length(indexOUT)){
      pond <-matrix(band_values,nrow = 1)%*%as.matrix(InRefl[indexOUT,])
      IntegrateChannel <-sum(band_values)
      Refl_Channel <-pond/IntegrateChannel
      OutRefl[[i]] <-matrix(Refl_Channel,nrow =1)
    } else {
      message('Please make sure that the spectral resolution for sensor response function ')
      message('__ is compatible ewith model (1 nm spectral sampling from 400 to 2500 m) __')
      message('______________________ This is currently not the case _____________________')

      message('_____ The follwing spectral band is not within accepted spectral range ____')
      print(SRF$Spectral_Bands[i])
      message('_________________________ Values will be set to 0  ________________________')
      OutRefl[[i]] <-matrix(0,ncol = ncol(InRefl),nrow =1)

    }
  }
  OutRefl <-do.call('rbind',OutRefl)
  return(OutRefl)
}

#' This function converts high resolution spectral info into broader spectral characteristics
#'
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for clear conditions
#' @param SRF list. Spectral response function, corresponding spectral bands, and Original Bands
#'
#' @return SpecSensor list. list of input specs with sensor resolution
#' @export
PrepareSensorSimulation <- function(SpecPROSPECT,SpecSOIL,SpecATM,SRF){

  # adjust optical constants
  wvl <- SpecPROSPECT$lambda
  # leaf properties
  SpecPROSPECT_Sensor = applySensorCharacteristics(wvl,SpecPROSPECT,SRF)
  SpecPROSPECT_Sensor = split(SpecPROSPECT_Sensor, rep(1:ncol(SpecPROSPECT_Sensor), each = nrow(SpecPROSPECT_Sensor)))
  names(SpecPROSPECT_Sensor)=names(SpecPROSPECT)
  # atmospheric properties
  SpecATM_Sensor = applySensorCharacteristics(wvl,SpecATM,SRF)
  SpecATM_Sensor = split(SpecATM_Sensor, rep(1:ncol(SpecATM_Sensor), each = nrow(SpecATM_Sensor)))
  names(SpecATM_Sensor)=names(SpecATM)
  # soil properties
  SpecSOIL_Sensor = applySensorCharacteristics(wvl,SpecSOIL,SRF)
  SpecSOIL_Sensor = split(SpecSOIL_Sensor, rep(1:ncol(SpecSOIL_Sensor), each = nrow(SpecSOIL_Sensor)))
  names(SpecSOIL_Sensor)=names(SpecSOIL)
  SpecSensor <- list('SpecPROSPECT_Sensor' = SpecPROSPECT_Sensor,'SpecATM_Sensor' = SpecATM_Sensor,'SpecSOIL_Sensor' = SpecSOIL_Sensor)
  return(list(SpecSensor))
}

#' Computes spectral response function based on wavelength and FWHM characteristics
#' @param wvl numeric. spectral sampling of the sensor
#' @param FWHM numeric. Full Width Half Maximum for each spectral band
#'
#' @return SRF list. Information about spectral response Spectral Bands of the sensor and Original Bands for which SRF is defined
#' @importFrom stats dnorm
#' @export
Compute_SRF <- function(wvl,FWHM){

  # define full spectral domain in optical domain
  lambda <- seq(400,2500,by = 1)
  VoidSpectrum <- matrix(0,nrow = length(lambda),ncol = 1)
  Spectral_Response <- matrix(0,ncol = length(wvl),nrow = length(lambda))
  for (i in 1:length(wvl)){
    y <- dnorm(lambda,wvl[i],FWHM[i]/2.355)
    Spectral_Response[,i] <- y/max(y)
  }
  Spectral_Response[which(Spectral_Response<0.001)] <- 0
  SRF <- list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=wvl,'OriginalBands'=lambda)
  return(SRF)
}

#' Computes spectral sensor characteristics for optical constants required for PROSPECT
#' soil parameters as well as direct and diffuse radiation for clear conditions
#' @param SensorName character. name of the sensor. Either already defined, or to be defined based on other inputs
#' @param SpectralProps list. Sensor spectral characteristics including wl = central wavelength and fwhm = corresponding FWHM
#' @param Path_SensorResponse character. Path where to get or save spectral response function
#'
#' @return SRF list. Spectral response function, corresponding spectral bands, and Original Bands
#' @export
get_spec_sensor <- function(SensorName = 'MyCustomSensor',
                            SpectralProps = NULL,
                            Path_SensorResponse = NULL){
  SRF <- GetRadiometry(SensorName,SpectralProps = SpectralProps,
                       Path_SensorResponse = Path_SensorResponse)
  if (is.null(SRF)){
    stop()
  }
  # adjust optical constants from 1nm sampling into spectral S2 spectral sampling
  SpecSensor <- PrepareSensorSimulation(prosail::SpecPROSPECT,
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
#'
#' @return List of S2 angles (SZA, SAA, VZA, VAA)
#' @importFrom XML xml xmlToList
#' @export
get_S2geometry <- function(MTD_TL_xml){
  # read XML file containing info about geometry of acquisition
  s2xml <- XML::xml(MTD_TL_xml)
  s2xml <- XML::xmlToList(s2xml)
  Distrib_SunAngle <- Distrib_ViewAngle <- list()
  # SZA
  Distrib_SunAngle$Zenith <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Zenith$Values_List
  SZA <- c()
  for (i in 1:length(Distrib_SunAngle$Zenith)){
    SZA <- c(SZA,as.numeric(strsplit(x = Distrib_SunAngle$Zenith[i]$VALUES,split = ' ')[[1]]))
  }
  SZA <- SZA[which(!is.na(SZA))]
  # SAA
  Distrib_SunAngle$Azimuth <- s2xml$Geometric_Info$Tile_Angles$Sun_Angles_Grid$Azimuth$Values_List
  SAA <- c()
  for (i in 1:length(Distrib_SunAngle$Zenith)){
    SAA <- c(SAA,as.numeric(strsplit(x = Distrib_SunAngle$Azimuth[i]$VALUES,split = ' ')[[1]]))
  }
  SAA <- SAA[which(!is.na(SAA))]
  Distrib_ViewAngle$Zenith <- s2xml$Geometric_Info$Tile_Angles$Viewing_Incidence_Angles_Grids$Zenith$Values_List
  VZA <- c()
  for (i in 1:length(Distrib_ViewAngle$Zenith)){
    VZA <- c(VZA,as.numeric(strsplit(x = Distrib_ViewAngle$Zenith[i]$VALUES,split = ' ')[[1]]))
  }
  VZA <- VZA[which(!is.na(VZA))]
  # VAA
  Distrib_ViewAngle$Azimuth <- s2xml$Geometric_Info$Tile_Angles$Viewing_Incidence_Angles_Grids$Azimuth$Values_List
  VAA <- c()
  for (i in 1:length(Distrib_ViewAngle$Azimuth)){
    VAA <- c(VAA,as.numeric(strsplit(x = Distrib_ViewAngle$Azimuth[i]$VALUES,split = ' ')[[1]]))
  }
  VAA <- VAA[which(!is.na(VAA))]
  return(list('SAA'=SAA,'SZA'=SZA,'VAA'=VAA,'VZA'=VZA))
}
