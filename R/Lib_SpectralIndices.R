# ==============================================================================
# prosail
# Lib_SpectralIndices.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ==============================================================================
# This Library includes aims at computing spectral indices from reflectance data
# ==============================================================================

#' this function aims at computing spectral indices from Sensor reflectance data in raster object
#' it computes the spectral indices based on their computation with Sentinel-2
#' and assumes that the bands of the S2 data follow this order
#' wavelength	= {496.6, 560.0, 664.5, 703.9, 740.2, 782.5, 835.1, 864.8, 1613.7, 2202.4}
#' Full description of the indices can be found here:
#' https://www.sentinel-hub.com/eotaxonomy/indices
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. wavelength in nanometers of the spectral bands of Refl.
#' @param Sel_Indices  list. list of spectral indices to be computed
#' @param StackOut logical. If TRUE returns a stack, otherwise a list of rasters.
#'
#' @return list. includes
#' - SpectralIndices: List of spectral indices computed from the reflectance initially provided
#' - listIndices: list of spectral indices computable with the function
#' @importFrom methods is
#' @import raster
#' @export

ComputeSpectralIndices_Raster <- function(Refl, SensorBands, Sel_Indices='ALL', StackOut=T){

  S2Bands <- c('B2'=496.6, 'B3'=560.0, 'B4'=664.5, 'B5'=703.9, 'B6'=740.2,
               'B7' = 782.5, 'B8' = 835.1, 'B8A' = 864.8, 'B11' = 1613.7, 'B12' = 2202.4)

  SpectralIndices <- list()
  Sen2S2 <- get_bands_close2s2(SensorBands,S2Bands)


  if(is.list(Refl))
    Refl <- raster::stack(Refl[Sen2S2]) # checks that all rasters have same crs/extent
  Refl <- subset(Refl, Sen2S2)
  if(is(Refl, 'RasterStack'))
    Refl = raster::brick(Refl) # loads the data: much faster than stack that is loading raster each time it is used
  if(!is(Refl, 'RasterBrick'))
    stop('Refl is expected to be a RasterStack or a list of rasters')

  names(Refl) <- names(Sen2S2)

  IndexAll <- list()
  # # set zero vaues to >0 in order to avoid problems
  # SelZero <- which(Refl==0)
  # Refl[SelZero] <- 0.005
  # if (dim(Refl)[1]==length(SensorBands)){
  #   Refl <- t(Refl)
  # }

  # inelegant but meeeeh
  listIndices <- list('ARI1','ARI2','ARVI','BAI','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDVI','NDVI_G',
                      'NDVI705','NDWI','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR')
  if (Sel_Indices[1]=='ALL'){
    Sel_Indices = listIndices
  }
  if ('ARI1'%in%Sel_Indices){
    ARI1 <- (1/Refl[["B3"]])-(1/Refl[["B5"]])
    SpectralIndices$ARI1 <- ARI1
  }
  if ('ARI2'%in%Sel_Indices){
    ARI2 <- (Refl[["B8"]]/Refl[["B2"]])-(Refl[["B8"]]/Refl[["B3"]])
    SpectralIndices$ARI2 <- ARI2
  }
  if ('ARVI'%in%Sel_Indices){
    ARVI <- (Refl[["B8"]]-(2*Refl[["B4"]]-Refl[["B2"]]))/(Refl[["B8"]]+(2*Refl[["B4"]]-Refl[["B2"]]))
    SpectralIndices$ARVI <- ARVI
  }
  if ('BAI'%in%Sel_Indices){
    BAI     = (1/((0.1-Refl[["B4"]])**2+(0.06-Refl[["B8"]])**2))
    SpectralIndices$BAI <- BAI
  }
  if ('CHL_RE'%in%Sel_Indices){
    CHL_RE  = Refl[["B5"]]/Refl[["B8"]]
    SpectralIndices$CHL_RE <- CHL_RE
  }
  if ('CRI1'%in%Sel_Indices){
    CRI1    = (1/Refl[["B2"]])-(1/Refl[["B3"]])
    SpectralIndices$CRI1 <- CRI1
  }
  if ('CRI2'%in%Sel_Indices){
    CRI2    = (1/Refl[["B2"]])-(1/Refl[["B5"]])
    SpectralIndices$CRI2 <- CRI2
  }
  if ('EVI'%in%Sel_Indices){
    EVI     = 2.5*(Refl[["B8"]]-Refl[["B4"]])/((Refl[["B8"]]+6*Refl[["B4"]]-(7.5*Refl[["B2"]]+1)))
    SpectralIndices$EVI <- EVI
  }
  if ('EVI2'%in%Sel_Indices){
    EVI2    = 2.5*(Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+2.4*Refl[["B4"]]+1)
    SpectralIndices$EVI2 <- EVI2
  }
  if ('GRVI1'%in%Sel_Indices){
    GRVI1   = (Refl[["B4"]]-Refl[["B3"]])/(Refl[["B4"]]+Refl[["B3"]])
    SpectralIndices$GRVI1 <- GRVI1
  }
  if ('GNDVI'%in%Sel_Indices){
    GNDVI   = (Refl[["B8"]]-Refl[["B3"]])/(Refl[["B8"]]+Refl[["B3"]])
    SpectralIndices$GNDVI <- GNDVI
  }
  if ('IRECI'%in%Sel_Indices){
    IRECI   = (Refl[["B7"]]-Refl[["B4"]])*(Refl[["B6"]]/(Refl[["B5"]]))
    SpectralIndices$IRECI <- IRECI
  }
  if ('LAI_SAVI'%in%Sel_Indices){
    LAI_SAVI = -log(0.371 + 1.5 * (Refl[["B8"]] - Refl[["B4"]]) / (Refl[["B8"]]+ Refl[["B4"]]+ 0.5)) / 2.4
    SpectralIndices$LAI_SAVI <- LAI_SAVI
  }
  if  ('MCARI'%in%Sel_Indices){
    MCARI   = (1-0.2*(Refl[["B5"]]-Refl[["B3"]])/(Refl[["B5"]]-Refl[["B4"]]))
    SpectralIndices$MCARI <- MCARI
  }
  if ('mNDVI705'%in%Sel_Indices){
    mNDVI705 = (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B6"]]+Refl[["B5"]]-2*Refl[["B2"]])
    SpectralIndices$mNDVI705 <- mNDVI705
  }
  if ('MSAVI2'%in%Sel_Indices){
    MSAVI2  = ((Refl[["B8"]]+1)-0.5*sqrt(((2*Refl[["B8"]])-1)**2+8*Refl[["B4"]]))
    SpectralIndices$MSAVI2 <- MSAVI2
  }
  if ('MSI'%in%Sel_Indices){
    MSI     = Refl[["B11"]]/Refl[["B8A"]]
    SpectralIndices$MSI <- MSI
  }
  if ('mSR705'%in%Sel_Indices){
    mSR705  = (Refl[["B6"]]-Refl[["B2"]])/(Refl[["B5"]]-Refl[["B2"]])
    SpectralIndices$mSR705 <- mSR705
  }
  if ('MTCI'%in%Sel_Indices){
    MTCI    = (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B5"]]+Refl[["B4"]])
    SpectralIndices$MTCI <- MTCI
  }
  if ('nBR_RAW'%in%Sel_Indices){
    nBR_RAW = (Refl[["B8"]]-Refl[["B12"]])/(Refl[["B8"]]+Refl[["B12"]])
    SpectralIndices$nBR_RAW <- nBR_RAW
  }
  if ('NDI_45'%in%Sel_Indices){
    NDI_45  = (Refl[["B5"]]-Refl[["B4"]])/(Refl[["B5"]]+Refl[["B4"]])
    SpectralIndices$NDI_45 <- NDI_45
  }
  if ('NDII'%in%Sel_Indices){
    NDII    = (Refl[["B8A"]]-Refl[["B11"]])/(Refl[["B8A"]]+Refl[["B11"]])
    SpectralIndices$NDII <- NDII
  }
  if ('NDVI'%in%Sel_Indices){
    NDVI   = (Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+Refl[["B4"]])
    SpectralIndices$NDVI <- NDVI
  }
  if ('NDVI_G'%in%Sel_Indices){
    NDVI_G  = Refl[["B3"]]*(Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+Refl[["B4"]])
    SpectralIndices$NDVI_G <- NDVI_G
  }
  if ('NDVI705'%in%Sel_Indices){
    NDVI705 = (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B6"]]+Refl[["B5"]])
    SpectralIndices$NDVI705 <- NDVI705
  }
  if ('NDWI'%in%Sel_Indices){
    NDWI    = (Refl[["B3"]]-Refl[["B8"]])/(Refl[["B3"]]+Refl[["B8"]])
    SpectralIndices$NDWI <- NDWI
  }
  if ('PSRI'%in%Sel_Indices){
    PSRI    = (Refl[["B4"]]-Refl[["B2"]])/(Refl[["B5"]])
    SpectralIndices$PSRI <- PSRI
  }
  if ('PSRI_NIR'%in%Sel_Indices){
    PSRI_NIR = (Refl[["B4"]]-Refl[["B2"]])/(Refl[["B8"]])
    SpectralIndices$PSRI_NIR <- PSRI_NIR
  }
  if ('RE_NDVI'%in%Sel_Indices){
    RE_NDVI = (Refl[["B8"]]-Refl[["B6"]])/(Refl[["B8"]]+Refl[["B6"]])
    SpectralIndices$RE_NDVI <- RE_NDVI
  }
  if ('RE_NDWI'%in%Sel_Indices){
    RE_NDWI = (Refl[["B4"]]-Refl[["B6"]])/(Refl[["B4"]]+Refl[["B6"]])
    SpectralIndices$RE_NDWI <- RE_NDWI
  }
  if ('S2REP'%in%Sel_Indices){
    S2REP   = 705+35*(0.5*(Refl[["B8"]]+Refl[["B5"]])-Refl[["B6"]])/(Refl[["B7"]]-Refl[["B6"]])
    SpectralIndices$S2REP <- S2REP
  }
  if ('SAVI'%in%Sel_Indices){
    SAVI    = 1.5*(Refl[["B8"]]-Refl[["B5"]])/(Refl[["B8"]]+Refl[["B5"]]+0.5)
    SpectralIndices$SAVI <- SAVI
  }
  if ('SIPI'%in%Sel_Indices){
    SIPI    = (Refl[["B8"]]-Refl[["B2"]])/(Refl[["B8"]]-Refl[["B4"]])
    SpectralIndices$SIPI <- SIPI
  }
  if ('SR'%in%Sel_Indices){
    SR <- Refl[["B8"]]/Refl[["B4"]]
    SpectralIndices$SR <- SR
  }
  if ('CR_SWIR'%in%Sel_Indices){
    CR_SWIR <- Refl[["B11"]]/(Refl[["B8A"]]+(S2Bands['B11']-S2Bands['B8A'])*(Refl[["B12"]]-Refl[["B8A"]])/(S2Bands['B12']-S2Bands['B8A']))
    SpectralIndices$CR_SWIR <- CR_SWIR
  }

  if(StackOut)
    SpectralIndices = raster::stack(SpectralIndices)

  res = list('SpectralIndices'=SpectralIndices,'listIndices'=listIndices)
  return(res)
}


#' this function identifies the bands of a given sensor with closest match to S2 bands
#'
#' @param SensorBands numeric. wavelength in nanometer of the sensor of interest
#' @param S2Bands numeric or list. Named vector or list of spectral bands corresponding to Sentinel-2
#'
#' @return numeric. band numbers of original sensor corresponding to S2
#' @export
get_bands_close2s2 <- function(SensorBands,S2Bands){
  sapply(S2Bands, function(x){b = which.min(abs(SensorBands-x)); names(b)=''; b})
}



# compute_lasrc_indices <- function(lasrc_dir){
#   # writes indices in files
#   require('raster')
#   source('R/functions.R')
#   LASRCBands <- c('band2'=496.6, 'band3'=560.0, 'band4'=664.5, 'band5'=703.9, 'band6'=740.2,
#                   'band7' = 782.5, 'band8' = 835.1, 'band11' = 1613.7, 'band12' = 2202.4, 'band8a' = 864.8)
#
#   bands_lasrc = data.frame(
#     name = c(sprintf('band%d', c(2:8, 11:12)), 'band8a'),
#     # res = c(rep(10, 3), rep(20, 3), 10, rep(20, 3)),
#     s2name = c(sprintf('B%02d', c(2:8, 11:12)), 'B08A'),
#     stringsAsFactors=F)
#
#   r = brick(stack(band_LASRC_file(lasrc_dir, bands_lasrc$name)))
#   indices = ComputeSpectralIndices_Raster(r, LASRCBands)
#   for(i in names(indices$SpectralIndices)){
#     file = file.path(lasrc_dir, paste0(basename(lasrc_dir), '_', i, '.tif'))
#     if(!file.exists(file)){
#       writeRaster(indices$SpectralIndices[[i]], filename = file)
#       cat(sprintf('File written: %s\n', file))
#     }
#   }
# }
