# ============================================================================== =
# prosail
# Lib_SpectralIndices.R
# ============================================================================== =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================== =
# This Library includes aims at computing spectral indices from reflectance data
# ============================================================================== =

#' This function computes Area under curve for continuum removed reflectances
#' See Malenovský et al. (2013) for details
#' http://dx.doi.org/10.1016/j.rse.2012.12.015
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param AUCminmax list. wavelengths of lower and upper boundaries ('CRmin' and 'CRmax')
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return AUCval raster
#' @export
AUC <- function(Refl, SensorBands, AUCminmax, ReflFactor=1){

  AUCbands <- list()
  AUCbands[['CRmin']] <- SensorBands[get_closest_bands(SensorBands,AUCminmax[['CRmin']])]
  AUCbands[['CRmax']] <- SensorBands[get_closest_bands(SensorBands,AUCminmax[['CRmax']])]
  Bands <- get_closest_bands(SensorBands,AUCbands)
  for (i in Bands[['CRmin']]:Bands[['CRmax']]){
    if (is.na(match(i,Bands))){
      AUCbands[[paste('B',i,sep = '')]] <- SensorBands[i]
    }
  }
  # compute continuum removal for all spectral bands
  CR <- CR_WL(Refl = Refl, SensorBands = SensorBands,
              CRbands = AUCbands, ReflFactor=ReflFactor)

  WL <- sort(unlist(AUCbands),decreasing = F)
  AUCval <- 0.5*(1-CR[[1]])*(WL[2]-WL[1])
  for (i in 2:length(CR)){
    AUCval <- AUCval+0.5*(2-CR[[i-1]]-CR[[i]])*(WL[i+1]-WL[i])
  }
  AUCval <- AUCval+0.5*(1-CR[[length(CR)]])*(WL[i+2]-WL[i+1])
  return(AUCval)
}

#' This function extracts boundaries to be used to compute continuum from reflectance data
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param CRbands list. list of spectral bands (central wavelength) including CRmin and CRmax
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return CRminmax list. list of rasters corresponding to minimum and maximum wavelengths
#' @export
CRbound <- function(Refl, SensorBands, CRbands, ReflFactor=1){

  # get closest spectral bands from CR1 and CR2
  Bands <- get_closest_bands(SensorBands,list(CRbands[['CRmin']],CRbands[['CRmax']]))
  WL <- SensorBands[Bands]
  # get equation for line going from CR1 to CR2
  CRminmax <- readRasterBands(Refl = Refl, Bands = Bands, ReflFactor=ReflFactor)
  names(CRminmax) <- paste('WL_',as.character(WL),sep = '')
  return(CRminmax)
}

#' This function extracts boundaries to be used to compute continuum from reflectance data
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param CRbands list. list of spectral bands (central wavelength) including CRmin and CRmax
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom progress progress_bar
#' @export
CR_WL <- function(Refl, SensorBands, CRbands, ReflFactor=1){

  # Make sure CRmin and CRmax are correctly defined
  if (is.na(match('CRmin',names(CRbands))) | is.na(match('CRmax',names(CRbands)))){
    stop('Please define CRmin and CRmax (CRmin<CRmax) as spectral bands in CRbands')
  }
  if (CRbands[['CRmax']] < CRbands[['CRmin']]){
    stop('Please define CRmin < CRmax in CRbands')
  }
  # extract CRmin and CRmax
  CRminmax <- CRbound(Refl, SensorBands, CRbands, ReflFactor=ReflFactor)
  # extract other bands and compute CR
  CRmin <- SensorBands[get_closest_bands(SensorBands,CRbands[['CRmin']])]
  CRmax <- SensorBands[get_closest_bands(SensorBands,CRbands[['CRmax']])]
  CRbands[['CRmin']] <- NULL
  CRbands[['CRmax']] <- NULL
  CR <- list()
  # initiate progress bar
  pgbarlength <- length(CRbands)
  pb <- progress_bar$new(
    format = "Computing continuum removal [:bar] :percent in :elapsedfull , estimated time remaining :eta",
    total = pgbarlength, clear = FALSE, width= 100)
  # computation for each band
  for (band in CRbands){
    pb$tick()
    bandrank <- get_closest_bands(SensorBands,band)
    raster2CR <- readRasterBands(Refl = Refl, Bands = bandrank, ReflFactor=ReflFactor)
    CR[[as.character(band)]] <- ComputeCR(WLmin = CRmin, WLmax = CRmax,
                                          WLtarget = band, boundaries=CRminmax,
                                          target=raster2CR)
    # CR[[as.character(band)]] <- target/(boundaries[[1]]+(WLtarget-WLmin)*(boundaries[[2]]-boundaries[[1]])/(WLmax-WLmin))
  }
  return(CR)
}

#' This function computes continuum removal value for a spectral band of interest,
#' based on lower and upper wavelengths corresponding to boundaries of the continuum
#'
#' @param WLmin numeric. wavelength of the spectral band corresponding to minimum boundary
#' @param WLmax numeric. wavelength of the spectral band corresponding to maximum boundary
#' @param WLtarget numeric. wavelength of the spectral band for which CR is computed
#' @param boundaries list. raster objects corresponding to minimum and maximum wavelengths
#' @param target list. raster object corresponding target wavelength
#'
#' @return CR list. raster object corresponding to continuum removed value
#' @export
ComputeCR <- function(WLmin, WLmax, WLtarget, boundaries, target){

  CR <- target/(boundaries[[1]]+(WLtarget-WLmin)*(boundaries[[2]]-boundaries[[1]])/(WLmax-WLmin))
  return(CR)
}

#' this function produces a spectral index from an expression defining a spectral index
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. wavelength in nanometers of the spectral bands of Refl.
#' @param ExpressIndex  character. expression corresponding to the spectral index to compute
#' @param listBands list. list of spectral bands defined in the 'ExpressIndex' variable
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param NameIndex character. name for the index to be computed, provided in the raster layer
#'
#' @return numeric. band numbers of original sensor corresponding to S2
#' @export
ComputeSpectralIndices_fromExpression <- function(Refl, SensorBands, ExpressIndex , listBands, ReflFactor=1,NameIndex = NULL){

  # define which bands to be used in the spectral index
  Bands <- get_closest_bands(SensorBands,listBands)

  ClassRaster <- class(Refl)
  if (ClassRaster=='RasterBrick' | ClassRaster=='RasterStack' | ClassRaster=='stars'){
    # if !ReflFactor == 1 then apply a reflectance factor
    if (ClassRaster=='stars'){
      Refl <- Refl[Bands]
    } else {
      Refl <- raster::subset(Refl, Bands)
    }
    if (!ReflFactor==1){
      Refl <- Refl/ReflFactor
    }
  } else if(is.list(Refl)){
    Refl <- raster::stack(Refl[Bands]) # checks that all rasters have same crs/extent
    if (!ReflFactor==1){
      Refl <- Refl/ReflFactor
    }
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object or a list of rasters')
  }
  names(Refl) <- gsub(pattern = 'B',replacement = 'Band',x = names(Bands))

  nbBands <- unique(as.numeric(gsub(pattern = 'B',
                                    replacement = '',
                                    x =  unlist(regmatches(ExpressIndex,
                                                           gregexpr("B[[:digit:]]+",
                                                                    ExpressIndex))))))
  sortBand <- sort(nbBands,index.return=T,decreasing = T)
  matches <- unique(unlist(regmatches(ExpressIndex, gregexpr("B[[:digit:]]+", ExpressIndex))))[sortBand$ix]
  replaces <- paste("Refl[['Band",gsub(pattern = 'B',replacement = '',x = matches),"']]",sep = '')
  ExpressIndex_Final <- ExpressIndex
  for (bb in 1:length(matches)){
    ExpressIndex_Final <- gsub(pattern = matches[bb], replacement = replaces[bb], x = ExpressIndex_Final)
  }
  SI <- eval(parse(text = ExpressIndex_Final))
  if (!is.null(NameIndex)){
    names(SI) <- NameIndex
  }
  return(SI)
}

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
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return list. includes
#' - SpectralIndices: List of spectral indices computed from the reflectance initially provided
#' - listIndices: list of spectral indices computable with the function
#' @importFrom methods is
#' @importFrom raster stack brick
#' @export

ComputeSpectralIndices_Raster <- function(Refl, SensorBands, Sel_Indices='ALL', StackOut=T,ReflFactor=1){

  S2Bands <- c('B2'=496.6, 'B3'=560.0, 'B4'=664.5, 'B5'=703.9, 'B6'=740.2,
               'B7' = 782.5, 'B8' = 835.1, 'B8A' = 864.8, 'B11' = 1613.7, 'B12' = 2202.4)

  SpectralIndices <- list()
  Sen2S2 <- get_closest_bands(SensorBands,S2Bands)
  ClassRaster <- class(Refl)
  if (ClassRaster=='RasterBrick' | ClassRaster=='RasterStack' | ClassRaster=='stars'){
    # if !ReflFactor == 1 then apply a reflectance factor
    if (ClassRaster=='stars'){
      Refl <- Refl[Sen2S2]
    } else {
      Refl <- raster::subset(Refl, Sen2S2)
    }
    if (!ReflFactor==1){
      Refl <- Refl/ReflFactor
    }
  } else if(is.list(Refl)){
    Refl <- raster::stack(Refl[Sen2S2]) # checks that all rasters have same crs/extent
    if (!ReflFactor==1){
      Refl <- Refl/ReflFactor
    }
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object or a list of rasters')
  }
  names(Refl) <- names(Sen2S2)

  IndexAll <- list()
  # # set zero vaues to >0 in order to avoid problems
  # SelZero <- which(Refl==0)
  # Refl[SelZero] <- 0.005
  # if (dim(Refl)[1]==length(SensorBands)){
  #   Refl <- t(Refl)
  # }

  # inelegant but meeeeh
  listIndices <- list('ARI1','ARI2','ARVI','BAI','BAIS2','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDVI','NDVI_G',
                      'NDVI705','NDWI1','NDWI2','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR')
  if (Sel_Indices[1]=='ALL'){
    Sel_Indices <- listIndices
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
    BAI <- (1/((0.1-Refl[["B4"]])**2+(0.06-Refl[["B8"]])**2))
    SpectralIndices$BAI <- BAI
  }
  if ('BAIS2'%in%Sel_Indices){
    BAIS2 <-  (1-((Refl[["B6"]]*Refl[["B7"]]*Refl[["B8A"]])/Refl[["B4"]])**0.5)*((Refl[["B12"]]-Refl[["B8A"]])/((Refl[["B12"]]+Refl[["B8A"]])**0.5)+1)
    SpectralIndices$BAIS2 <- BAIS2
  }
  if ('CHL_RE'%in%Sel_Indices){
    CHL_RE <- Refl[["B5"]]/Refl[["B8"]]
    SpectralIndices$CHL_RE <- CHL_RE
  }
  if ('CRI1'%in%Sel_Indices){
    CRI1 <- (1/Refl[["B2"]])-(1/Refl[["B3"]])
    SpectralIndices$CRI1 <- CRI1
  }
  if ('CRI2'%in%Sel_Indices){
    CRI2 <- (1/Refl[["B2"]])-(1/Refl[["B5"]])
    SpectralIndices$CRI2 <- CRI2
  }
  if ('EVI'%in%Sel_Indices){
    EVI <- 2.5*(Refl[["B8"]]-Refl[["B4"]])/((Refl[["B8"]]+6*Refl[["B4"]]-7.5*Refl[["B2"]]+1))
    SpectralIndices$EVI <- EVI
  }
  if ('EVI2'%in%Sel_Indices){
    EVI2 <- 2.5*(Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+2.4*Refl[["B4"]]+1)
    SpectralIndices$EVI2 <- EVI2
  }
  if ('GRVI1'%in%Sel_Indices){
    GRVI1 <- (Refl[["B4"]]-Refl[["B3"]])/(Refl[["B4"]]+Refl[["B3"]])
    SpectralIndices$GRVI1 <- GRVI1
  }
  if ('GNDVI'%in%Sel_Indices){
    GNDVI <- (Refl[["B8"]]-Refl[["B3"]])/(Refl[["B8"]]+Refl[["B3"]])
    SpectralIndices$GNDVI <- GNDVI
  }
  if ('IRECI'%in%Sel_Indices){
    IRECI <- (Refl[["B7"]]-Refl[["B4"]])*(Refl[["B6"]]/(Refl[["B5"]]))
    SpectralIndices$IRECI <- IRECI
  }
  if ('LAI_SAVI'%in%Sel_Indices){
    LAI_SAVI <- -log(0.371 + 1.5 * (Refl[["B8"]] - Refl[["B4"]]) / (Refl[["B8"]]+ Refl[["B4"]]+ 0.5)) / 2.4
    SpectralIndices$LAI_SAVI <- LAI_SAVI
  }
  if  ('MCARI'%in%Sel_Indices){
    MCARI <- (1-0.2*(Refl[["B5"]]-Refl[["B3"]])/(Refl[["B5"]]-Refl[["B4"]]))
    SpectralIndices$MCARI <- MCARI
  }
  if ('mNDVI705'%in%Sel_Indices){
    mNDVI705 <- (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B6"]]+Refl[["B5"]]-2*Refl[["B2"]])
    SpectralIndices$mNDVI705 <- mNDVI705
  }
  if ('MSAVI2'%in%Sel_Indices){
    MSAVI2 <- ((Refl[["B8"]]+1)-0.5*sqrt(((2*Refl[["B8"]])-1)**2+8*Refl[["B4"]]))
    SpectralIndices$MSAVI2 <- MSAVI2
  }
  if ('MSI'%in%Sel_Indices){
    MSI <- Refl[["B11"]]/Refl[["B8A"]]
    SpectralIndices$MSI <- MSI
  }
  if ('mSR705'%in%Sel_Indices){
    mSR705 <- (Refl[["B6"]]-Refl[["B2"]])/(Refl[["B5"]]-Refl[["B2"]])
    SpectralIndices$mSR705 <- mSR705
  }
  if ('MTCI'%in%Sel_Indices){
    MTCI <- (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B5"]]+Refl[["B4"]])
    SpectralIndices$MTCI <- MTCI
  }
  if ('nBR_RAW'%in%Sel_Indices){
    nBR_RAW <- (Refl[["B8"]]-Refl[["B12"]])/(Refl[["B8"]]+Refl[["B12"]])
    SpectralIndices$nBR_RAW <- nBR_RAW
  }
  if ('NDI_45'%in%Sel_Indices){
    NDI_45 <- (Refl[["B5"]]-Refl[["B4"]])/(Refl[["B5"]]+Refl[["B4"]])
    SpectralIndices$NDI_45 <- NDI_45
  }
  if ('NDII'%in%Sel_Indices){
    NDII <- (Refl[["B8A"]]-Refl[["B11"]])/(Refl[["B8A"]]+Refl[["B11"]])
    SpectralIndices$NDII <- NDII
  }
  if ('NDVI'%in%Sel_Indices){
    NDVI <- (Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+Refl[["B4"]])
    SpectralIndices$NDVI <- NDVI
  }
  if ('NDVI_G'%in%Sel_Indices){
    NDVI_G <- Refl[["B3"]]*(Refl[["B8"]]-Refl[["B4"]])/(Refl[["B8"]]+Refl[["B4"]])
    SpectralIndices$NDVI_G <- NDVI_G
  }
  if ('NDVI705'%in%Sel_Indices){
    NDVI705 <- (Refl[["B6"]]-Refl[["B5"]])/(Refl[["B6"]]+Refl[["B5"]])
    SpectralIndices$NDVI705 <- NDVI705
  }
  if ('NDWI1'%in%Sel_Indices){
    NDWI1 <- (Refl[["B8A"]]-Refl[["B11"]])/(Refl[["B8A"]]+Refl[["B11"]])
    SpectralIndices$NDWI1 <- NDWI1
  }
  if ('NDWI2'%in%Sel_Indices){
    NDWI2 <- (Refl[["B8A"]]-Refl[["B12"]])/(Refl[["B8A"]]+Refl[["B12"]])
    SpectralIndices$NDWI2 <- NDWI2
  }
  if ('PSRI'%in%Sel_Indices){
    PSRI <- (Refl[["B4"]]-Refl[["B2"]])/(Refl[["B5"]])
    SpectralIndices$PSRI <- PSRI
  }
  if ('PSRI_NIR'%in%Sel_Indices){
    PSRI_NIR <- (Refl[["B4"]]-Refl[["B2"]])/(Refl[["B8"]])
    SpectralIndices$PSRI_NIR <- PSRI_NIR
  }
  if ('RE_NDVI'%in%Sel_Indices){
    RE_NDVI <- (Refl[["B8"]]-Refl[["B6"]])/(Refl[["B8"]]+Refl[["B6"]])
    SpectralIndices$RE_NDVI <- RE_NDVI
  }
  if ('RE_NDWI'%in%Sel_Indices){
    RE_NDWI <- (Refl[["B4"]]-Refl[["B6"]])/(Refl[["B4"]]+Refl[["B6"]])
    SpectralIndices$RE_NDWI <- RE_NDWI
  }
  if ('S2REP'%in%Sel_Indices){
    S2REP <- 705+35*(0.5*(Refl[["B8"]]+Refl[["B5"]])-Refl[["B6"]])/(Refl[["B7"]]-Refl[["B6"]])
    SpectralIndices$S2REP <- S2REP
  }
  if ('SAVI'%in%Sel_Indices){
    SAVI <- 1.5*(Refl[["B8"]]-Refl[["B5"]])/(Refl[["B8"]]+Refl[["B5"]]+0.5)
    SpectralIndices$SAVI <- SAVI
  }
  if ('SIPI'%in%Sel_Indices){
    SIPI <- (Refl[["B8"]]-Refl[["B2"]])/(Refl[["B8"]]-Refl[["B4"]])
    SpectralIndices$SIPI <- SIPI
  }
  if ('SR'%in%Sel_Indices){
    SR <- Refl[["B8"]]/Refl[["B4"]]
    SpectralIndices$SR <- SR
  }
  if ('TCARI'%in%Sel_Indices){
    SR <-

      Refl[["B8"]]/Refl[["B4"]]
    SpectralIndices$SR <- SR
  }
  if ('CR_SWIR'%in%Sel_Indices){
    CR_SWIR <- Refl[["B11"]]/(Refl[["B8A"]]+(S2Bands['B11']-S2Bands['B8A'])*(Refl[["B12"]]-Refl[["B8A"]])/(S2Bands['B12']-S2Bands['B8A']))
    SpectralIndices$CR_SWIR <- CR_SWIR
  }

  if(StackOut)
    SpectralIndices <- raster::stack(SpectralIndices)

  res <- list('SpectralIndices'=SpectralIndices,'listIndices'=listIndices)
  return(res)
}

#' this function aims at computing spectral indices from Sensor reflectance data.
#' it computes the spectral indices based on their computation with Sentinel-2
#' and assumes that the bands of the S2 data follow this order
#' wavelength	= {496.6, 560.0, 664.5, 703.9, 740.2, 782.5, 835.1, 864.8, 1613.7, 2202.4}
#' Full description of the indices can be found here:
#' https://www.sentinel-hub.com/eotaxonomy/indices
#'
#' @param Refl numeric. Reflectance dataset defined in matrix
#' @param Sel_Indices list. list of spectral indices to be computed
#' @param SensorBands numeric. wavelength of the spectral bands corresponding to the spectral index
#'
#' @return list. includes
#' - SpectralIndices: List of spectral indices computed from the reflectance initially provided
#' - listIndices: list of spectral indices computable with the function
#' @export

ComputeSpectralIndices_HS <- function(Refl,SensorBands,Sel_Indices='ALL'){

  SpectralIndices <- list()
  S2Bands <- data.frame('B2'=496.6, 'B3'=560.0, 'B4'=664.5, 'B5'=703.9, 'B6'=740.2,
                        'B7' = 782.5, 'B8' = 835.1, 'B8A' = 864.8, 'B11' = 1613.7, 'B12' = 2202.4)
  Sen2S2 <- get_closest_bands(SensorBands,S2Bands)
  IndexAll <- list()
  # set zero vaues to >0 in order to avoid problems
  SelZero <- which(Refl==0)
  Refl[SelZero] <- 0.005
  if (dim(Refl)[1]==length(SensorBands)){
    Refl <- t(Refl)
  }

  # inelegant but meeeeh
  listIndices <- list('ARI1','ARI2','ARVI','BAI','BAIS2','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDVI','NDVI_G',
                      'NDVI705','NDWI1','NDWI2','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR')
  if (Sel_Indices=='ALL'){
    Sel_Indices = listIndices
  }
  if ('ARI1'%in%Sel_Indices){
    ARI1 <- (1/Refl[,Sen2S2[["B3"]]])-(1/Refl[,Sen2S2[["B5"]]])
    SpectralIndices$ARI1 <- ARI1
  }
  if ('ARI2'%in%Sel_Indices){
    ARI2 <- (Refl[,Sen2S2[["B8"]]]/Refl[,Sen2S2[["B2"]]])-(Refl[,Sen2S2[["B8"]]]/Refl[,Sen2S2[["B3"]]])
    SpectralIndices$ARI2 <- ARI2
  }
  if ('ARVI'%in%Sel_Indices){
    ARVI <- (Refl[,Sen2S2[["B8"]]]-(2*Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B2"]]]))/(Refl[,Sen2S2[["B8"]]]+(2*Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B2"]]]))
    SpectralIndices$ARVI <- ARVI
  }
  if ('BAI'%in%Sel_Indices){
    BAI <- (1/((0.1-Refl[,Sen2S2[["B4"]]])**2+(0.06-Refl[,Sen2S2[["B8"]]])**2))
    SpectralIndices$BAI <- BAI
  }
  if ('BAIS2'%in%Sel_Indices){
    BAIS2 <-  (1-((Refl[,Sen2S2[["B6"]]]*Refl[,Sen2S2[["B7"]]]*Refl[,Sen2S2[["B8A"]]])/Refl[,Sen2S2[["B4"]]])**0.5)*((Refl[,Sen2S2[["B12"]]]-Refl[,Sen2S2[["B8A"]]])/((Refl[,Sen2S2[["B12"]]]+Refl[,Sen2S2[["B8A"]]])**0.5)+1)
    SpectralIndices$BAIS2 <- BAIS2
  }
  if ('CHL_RE'%in%Sel_Indices){
    CHL_RE <- Refl[,Sen2S2[["B5"]]]/Refl[,Sen2S2[["B8"]]]
    SpectralIndices$CHL_RE <- CHL_RE
  }
  if ('CRI1'%in%Sel_Indices){
    CRI1 <- (1/Refl[,Sen2S2[["B2"]]])-(1/Refl[,Sen2S2[["B3"]]])
    SpectralIndices$CRI1 <- CRI1
  }
  if ('CRI2'%in%Sel_Indices){
    CRI2 <- (1/Refl[,Sen2S2[["B2"]]])-(1/Refl[,Sen2S2[["B5"]]])
    SpectralIndices$CRI2 <- CRI2
  }
  if ('EVI'%in%Sel_Indices){
    EVI <- 2.5*(Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B4"]]])/((Refl[,Sen2S2[["B8"]]]+6*Refl[,Sen2S2[["B4"]]]-7.5*Refl[,Sen2S2[["B2"]]]+1))
    SpectralIndices$EVI <- EVI
  }
  if ('EVI2'%in%Sel_Indices){
    EVI2 <- 2.5*(Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B4"]]])/(Refl[,Sen2S2[["B8"]]]+2.4*Refl[,Sen2S2[["B4"]]]+1)
    SpectralIndices$EVI2 <- EVI2
  }
  if ('GRVI1'%in%Sel_Indices){
    GRVI1 <- (Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B3"]]])/(Refl[,Sen2S2[["B4"]]]+Refl[,Sen2S2[["B3"]]])
    SpectralIndices$GRVI1 <- GRVI1
  }
  if ('GNDVI'%in%Sel_Indices){
    GNDVI <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B3"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B3"]]])
    SpectralIndices$GNDVI <- GNDVI
  }
  if ('IRECI'%in%Sel_Indices){
    IRECI <- (Refl[,Sen2S2[["B7"]]]-Refl[,Sen2S2[["B4"]]])*(Refl[,Sen2S2[["B6"]]]/(Refl[,Sen2S2[["B5"]]]))
    SpectralIndices$IRECI <- IRECI
  }
  if ('LAI_SAVI'%in%Sel_Indices){
    LAI_SAVI <- -log(0.371 + 1.5 * (Refl[,Sen2S2[["B8"]]] - Refl[,Sen2S2[["B4"]]]) / (Refl[,Sen2S2[["B8"]]]+ Refl[,Sen2S2[["B4"]]]+ 0.5)) / 2.4
    SpectralIndices$LAI_SAVI <- LAI_SAVI
  }
  if  ('MCARI'%in%Sel_Indices){
    MCARI <- (1-0.2*(Refl[,Sen2S2[["B5"]]]-Refl[,Sen2S2[["B3"]]])/(Refl[,Sen2S2[["B5"]]]-Refl[,Sen2S2[["B4"]]]))
    SpectralIndices$MCARI <- MCARI
  }
  if ('mNDVI705'%in%Sel_Indices){
    mNDVI705 <- (Refl[,Sen2S2[["B6"]]]-Refl[,Sen2S2[["B5"]]])/(Refl[,Sen2S2[["B6"]]]+Refl[,Sen2S2[["B5"]]]-2*Refl[,Sen2S2[["B2"]]])
    SpectralIndices$mNDVI705 <- mNDVI705
  }
  if ('MSAVI2'%in%Sel_Indices){
    MSAVI2 <- ((Refl[,Sen2S2[["B8"]]]+1)-0.5*sqrt(((2*Refl[,Sen2S2[["B8"]]])-1)**2+8*Refl[,Sen2S2[["B4"]]]))
    SpectralIndices$MSAVI2 <- MSAVI2
  }
  if ('MSI'%in%Sel_Indices){
    MSI <- Refl[,Sen2S2[["B11"]]]/Refl[,Sen2S2[["B8"]]]
    SpectralIndices$MSI <- MSI
  }
  if ('mSR705'%in%Sel_Indices){
    mSR705 <- (Refl[,Sen2S2[["B6"]]]-Refl[,Sen2S2[["B2"]]])/(Refl[,Sen2S2[["B5"]]]-Refl[,Sen2S2[["B2"]]])
    SpectralIndices$mSR705 <- mSR705
  }
  if ('MTCI'%in%Sel_Indices){
    MTCI <- (Refl[,Sen2S2[["B6"]]]-Refl[,Sen2S2[["B5"]]])/(Refl[,Sen2S2[["B5"]]]+Refl[,Sen2S2[["B4"]]])
    SpectralIndices$MTCI <- MTCI
  }
  if ('nBR_RAW'%in%Sel_Indices){
    nBR_RAW <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B12"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B12"]]])
    SpectralIndices$nBR_RAW <- nBR_RAW
  }
  if ('NDI_45'%in%Sel_Indices){
    NDI_45 <- (Refl[,Sen2S2[["B5"]]]-Refl[,Sen2S2[["B4"]]])/(Refl[,Sen2S2[["B5"]]]+Refl[,Sen2S2[["B4"]]])
    SpectralIndices$NDI_45 <- NDI_45
  }
  if ('NDII'%in%Sel_Indices){
    NDII <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B11"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B11"]]])
    SpectralIndices$NDII <- NDII
  }
  if ('NDVI'%in%Sel_Indices){
    NDVI <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B4"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B4"]]])
    SpectralIndices$NDVI <- NDVI
  }
  if ('NDVI_G'%in%Sel_Indices){
    NDVI_G <- Refl[,Sen2S2[["B3"]]]*(Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B4"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B4"]]])
    SpectralIndices$NDVI_G <- NDVI_G
  }
  if ('NDVI705'%in%Sel_Indices){
    NDVI705 <- (Refl[,Sen2S2[["B6"]]]-Refl[,Sen2S2[["B5"]]])/(Refl[,Sen2S2[["B6"]]]+Refl[,Sen2S2[["B5"]]])
    SpectralIndices$NDVI705 <- NDVI705
  }
  if ('NDWI1'%in%Sel_Indices){
    NDWI1 <- (Refl[,Sen2S2[["B8A"]]]-Refl[,Sen2S2[["B11"]]])/(Refl[,Sen2S2[["B8A"]]]+Refl[,Sen2S2[["B11"]]])
    SpectralIndices$NDWI1 <- NDWI1
  }
  if ('NDWI2'%in%Sel_Indices){
    NDWI2 <- (Refl[,Sen2S2[["B8A"]]]-Refl[,Sen2S2[["B12"]]])/(Refl[,Sen2S2[["B8A"]]]+Refl[,Sen2S2[["B12"]]])
    SpectralIndices$NDWI2 <- NDWI2
  }
  if ('PSRI'%in%Sel_Indices){
    PSRI <- (Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B2"]]])/(Refl[,Sen2S2[["B5"]]])
    SpectralIndices$PSRI <- PSRI
  }
  if ('PSRI_NIR'%in%Sel_Indices){
    PSRI_NIR <- (Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B2"]]])/(Refl[,Sen2S2[["B8"]]])
    SpectralIndices$PSRI_NIR <- PSRI_NIR
  }
  if ('RE_NDVI'%in%Sel_Indices){
    RE_NDVI <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B6"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B6"]]])
    SpectralIndices$RE_NDVI <- RE_NDVI
  }
  if ('RE_NDWI'%in%Sel_Indices){
    RE_NDWI <- (Refl[,Sen2S2[["B4"]]]-Refl[,Sen2S2[["B6"]]])/(Refl[,Sen2S2[["B4"]]]+Refl[,Sen2S2[["B6"]]])
    SpectralIndices$RE_NDWI <- RE_NDWI
  }
  if ('S2REP'%in%Sel_Indices){
    S2REP <- 705+35*(0.5*(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B5"]]])-Refl[,Sen2S2[["B6"]]])/(Refl[,Sen2S2[["B7"]]]-Refl[,Sen2S2[["B6"]]])
    SpectralIndices$S2REP <- S2REP
  }
  if ('SAVI'%in%Sel_Indices){
    SAVI <- 1.5*(Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B5"]]])/(Refl[,Sen2S2[["B8"]]]+Refl[,Sen2S2[["B5"]]]+0.5)
    SpectralIndices$SAVI <- SAVI
  }
  if ('SIPI'%in%Sel_Indices){
    SIPI <- (Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B2"]]])/(Refl[,Sen2S2[["B8"]]]-Refl[,Sen2S2[["B4"]]])
    SpectralIndices$SIPI <- SIPI
  }
  if ('SR'%in%Sel_Indices){
    SR <- Refl[,Sen2S2[["B8"]]]/Refl[,Sen2S2[["B4"]]]
    SpectralIndices$SR <- SR
  }
  if ('CR_SWIR'%in%Sel_Indices){
    CR_SWIR <- Refl[,Sen2S2[["B11"]]]/(Refl[,Sen2S2[["B8A"]]]+(S2Bands$B11-S2Bands$B8A)*(Refl[,Sen2S2[["B12"]]]-Refl[,Sen2S2[["B8A"]]])/(S2Bands$B12-S2Bands$B8A))
    SpectralIndices$CR_SWIR <- CR_SWIR
  }
  res <- list('SpectralIndices'=SpectralIndices,'listIndices'=listIndices)
  return(res)
}

#' this function identifies the bands of a given sensor with closest match to its spectral characteristics
#'
#' @param SensorBands numeric. wavelength in nanometer of the sensor of interest
#' @param listBands numeric or list. Named vector or list of spectral bands corresponding to sensor
#'
#' @return numeric. band numbers of original sensor
#' @export
get_closest_bands <- function(SensorBands,listBands){
  sapply(listBands, function(x){b = which.min(abs(SensorBands-x)); names(b)=''; b})
}

#' This function computes interquartile range (IQR) criterion, which can be used
#' as a criterion for outlier detection
#'
#' @param DistVal numeric. vector of distribution of values
#' @param weightIRQ numeric. weighting factor appplied to IRQ to define lower and upper boudaries for outliers
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom stats IQR quantile
#' @export
IQR_outliers <- function(DistVal,weightIRQ = 1.5){
  iqr <- IQR(DistVal, na.rm=TRUE)
  range_IQR <- c(quantile(DistVal, 0.25,na.rm=TRUE),quantile(DistVal, 0.75,na.rm=TRUE))
  outlier_IQR <- c(range_IQR[1]-weightIRQ*iqr,range_IQR[2]+weightIRQ*iqr)
  return(outlier_IQR)
}

#' This function selects bands from a raster or stars object
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param Bands numeric. rank of bands to be read in Refl
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return Robj list. R object (default is raster, stars if Refl is stars object)
#' @importFrom raster subset stack
#' @export
readRasterBands <- function(Refl, Bands, ReflFactor=1){

  # get equation for line going from CR1 to CR2
  ClassRaster <- class(Refl)
  if (ClassRaster=='RasterBrick' | ClassRaster=='RasterStack' | ClassRaster=='stars'){
    # if !ReflFactor == 1 then apply a reflectance factor
    if (ClassRaster=='stars'){
      Robj <- Refl[Bands]
    } else {
      Robj <- raster::subset(Refl, Bands)
    }
    if (!ReflFactor==1){
      Robj <- Robj/ReflFactor
    }
  } else if(is.list(Refl)){
    Robj <- raster::stack(Refl[Bands]) # checks that all rasters have same crs/extent
    if (!ReflFactor==1){
      Robj <- Robj/ReflFactor
    }
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object or a list of rasters')
  }
  return(Robj)
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
