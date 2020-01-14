# ==============================================================================
# prosail
# Lib_SENSOR.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2019/11 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to convrsion from high resolution
# reflectance to sensor reflectance
# using either known snsor response, or gaussian filters
# ==============================================================================

#' reads spectral response from known sensor
#' spectral response from Sentinel-2 is already defined
#' @param SensorName character. name of the sensor
#' @param Path_SensorResponse character.
#' @return SRF list. Spectral response function, corresponding spectral bands, and Original Bands
#' @importFrom utils read.csv
#'
#' @export
GetRadiometry <- function(SensorName,Path_SensorResponse = NULL){

  # if sensor is SENTINEL-2
  if (SensorName=='Sentinel_2'){
    Spectral_Response <- Sentinel_2$Spectral_Response
    Spectral_Bands <- Sentinel_2$Spectral_Bands
    OriginalBands <- Sentinel_2$OriginalBands
  } else if (SensorName=='Sentinel_2A'){
    Spectral_Response <- Sentinel_2A$Spectral_Response
    Spectral_Bands <- Sentinel_2A$Spectral_Bands
    OriginalBands <- Sentinel_2A$OriginalBands
  } else if (SensorName=='Sentinel_2B'){
    Spectral_Response <- Sentinel_2B$Spectral_Response
    Spectral_Bands <- Sentinel_2B$Spectral_Bands
    OriginalBands <- Sentinel_2B$OriginalBands
    # if sensor is not SENTINEL-2
  }  else {
    # identify file containing spectral response
    if (!is.null(Path_SensorResponse)){
      Path_SRF <- file.path(Path_SensorResponse,paste(SensorName,'_Spectral_Response.csv',sep=''))
    } else {
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
    } else {
      message('___ Spectral response of the sensor expected here: ____')
      print(Path_SRF)
      message('_ Please check how to define sensor spectral response _')
      message('_____________ and create file accordingly _____________')
      stop()
    }
    Spectral_Response <- t(Spectral_Response)
  }
  SRF = list("Spectral_Response"=Spectral_Response, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)
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
      stop()
    }
  }
  OutRefl <-do.call('rbind',OutRefl)
  return(OutRefl)
}
