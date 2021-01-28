# ============================================================================= =
# prosail
# Lib_PROSAIL_HybridInversion.R
# ============================================================================= =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================= =
# This Library includes functions dedicated to PROSAIL inversion using hybrid
# approach based on SVM regression
# ============================================================================= =

#' This function performs full training for hybrid invrsion using SVR with
#' values for default parameters
#'
#' @param minval list. minimum value for input parameters sampled to produce a training LUT
#' @param maxval list. maximum value for input parameters sampled to produce a training LUT
#' @param TypeDistrib  list. Type of distribution. Either 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean value and STD corresponding to the parameters sampled with gaussian distribution
#' @param ParmSet list. list of input parameters set to a specific value
#' @param nbSamples numeric. number of samples in training LUT
#' @param nbSamplesPerRun numeric. number of training sample per individual regression model
#' @param nbModels numeric. number of individual models to be run for ensemble
#' @param Replacement bolean. is there replacement in subsampling?
#' @param SAILversion character. Either 4SAIL or 4SAIL2
#' @param Parms2Estimate list. list of input parameters to be estimated
#' @param Bands2Select list. list of bands used for regression for each input parameter
#' @param NoiseLevel list. list of noise value added to reflectance (defined per input parm)
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for clear conditions
#' @param Path_Results character. path for results
#'
#' @return modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
#' @export

train_prosail_inversion <- function(minval=NULL,maxval=NULL,
                                    TypeDistrib=NULL,GaussianDistrib= NULL,ParmSet=NULL,
                                    nbSamples=2000,nbSamplesPerRun=100,nbModels=20,Replacement=TRUE,
                                    SAILversion='4SAIL',
                                    Parms2Estimate='lai',Bands2Select=NULL,NoiseLevel=NULL,
                                    SpecPROSPECT = NULL, SpecSOIL = NULL, SpecATM = NULL,
                                    Path_Results='./'){

  #########################################################################
  ###           1- PRODUCE A LUT TO TRAIN THE HYBRID INVERSION          ###
  #########################################################################
  # Define sensor characteristics
  if (is.null(SpecPROSPECT)){
    SpecPROSPECT <- prosail::SpecPROSPECT
  }
  if (is.null(SpecSOIL)){
    SpecSOIL <- prosail::SpecSOIL
  }
  if (is.null(SpecPROSPECT)){
    SpecATM <- prosail::SpecATM
  }
  # define distribution for parameters to be sampled
  if (is.null(TypeDistrib)){
    TypeDistrib <- data.frame('CHL'='Uniform', 'CAR'='Uniform','EWT' = 'Uniform','ANT' = 'Uniform','LMA' = 'Uniform','N' = 'Uniform', 'BROWN'='Uniform',
                              'psoil' = 'Uniform','LIDFa' = 'Uniform', 'lai' = 'Uniform','q'='Uniform','tto' = 'Uniform','tts' = 'Uniform', 'psi' = 'Uniform')
  }
  if (is.null(GaussianDistrib)){
    GaussianDistrib <- list('Mean'=NULL,'Std'=NULL)
  }
  if (is.null(minval)){
    minval <- data.frame('CHL'=10,'CAR'=0,'EWT' = 0.01,'ANT' = 0,'LMA' = 0.005,'N' = 1.0,'psoil' = 0.0, 'BROWN'=0.0,
                         'LIDFa' = 20, 'lai' = 0.5,'q'=0.1,'tto' = 0,'tts' = 20, 'psi' = 80)
  }
  if (is.null(maxval)){
    maxval <- data.frame('CHL'=75,'CAR'=15,'EWT' = 0.03,'ANT' = 2,'LMA' = 0.03,'N' = 2.0, 'psoil' = 1.0, 'BROWN'=0.5,
                         'LIDFa' = 70, 'lai' = 7,'q'=0.2,'tto' = 5,'tts' = 30, 'psi' = 110)
  }
  # define min and max values
  # fixed parameters
  if (is.null(ParmSet)){
    ParmSet <- data.frame('TypeLidf' = 2, 'alpha' = 40)
  }
  # produce input parameters distribution
  InputPROSAIL <- get_distribution_input_prosail(minval,maxval,ParmSet,nbSamples,
                                                 TypeDistrib = TypeDistrib,
                                                 Mean = GaussianDistrib$Mean,Std = GaussianDistrib$Std)
  if (SAILversion=='4SAIL2'){
    # Definition of Cv & update LAI
    MaxLAI <- min(c(maxval$lai),4)
    InputPROSAIL$Cv <- NA*InputPROSAIL$lai
    InputPROSAIL$Cv[which(InputPROSAIL$lai>MaxLAI)] <- 1
    InputPROSAIL$Cv[which(InputPROSAIL$lai<=MaxLAI)] <- (1/MaxLAI)+InputPROSAIL$lai[which(InputPROSAIL$lai<=MaxLAI)]/(MaxLAI+1)
    InputPROSAIL$Cv <- InputPROSAIL$Cv*matrix(rnorm(length(InputPROSAIL$Cv),mean = 1,sd = 0.1))
    InputPROSAIL$Cv[which(InputPROSAIL$Cv<0)] <- 0
    InputPROSAIL$Cv[which(InputPROSAIL$Cv>1)] <- 1
    InputPROSAIL$Cv[which(InputPROSAIL$lai>MaxLAI)] <- 1
    InputPROSAIL$fraction_brown <- 0+0*InputPROSAIL$lai
    InputPROSAIL$diss <- 0+0*InputPROSAIL$lai
    InputPROSAIL$Zeta <- 0.2+0*InputPROSAIL$lai
    InputPROSAIL$lai <- InputPROSAIL$lai*InputPROSAIL$Cv
  }

  # generate LUT of BRF corresponding to InputPROSAIL, for a sensor
  BRF_LUT <- Generate_LUT_BRF(SAILversion=SAILversion,InputPROSAIL = InputPROSAIL,
                              SpecPROSPECT = SpecPROSPECT,SpecSOIL = SpecSOIL,SpecATM = SpecATM)

  # write parameters LUT
  output <- matrix(unlist(InputPROSAIL), ncol = length(InputPROSAIL), byrow = FALSE)
  filename <- file.path(Path_Results,'PROSAIL_LUT_InputParms.txt')
  write.table(x = format(output, digits=3),file = filename,append = F, quote = F,
              col.names = names(InputPROSAIL), row.names = F,sep = '\t')
  # Write BRF LUT corresponding to parameters LUT
  filename <- file.path(Path_Results,'PROSAIL_LUT_Reflectance.txt')
  write.table(x = format(t(BRF_LUT), digits=5),file = filename,append = F, quote = F,
              col.names = SpecPROSPECT$lambda, row.names = F,sep = '\t')

  # Which bands will be used for inversion?
  if (is.null(Bands2Select)){
    Bands2Select <- list()
    for (parm in Parms2Estimate){
      Bands2Select[[parm]] <- seq(1,length(SpecPROSPECT$lambda))
    }
  }
  # Add gaussian noise to reflectance LUT: one specific LUT per parameter
  if (is.null(NoiseLevel)){
    NoiseLevel <- list()
    for (parm in Parms2Estimate){
      NoiseLevel[[parm]] <- 0.01
    }
  }

  # produce LIT with noise
  BRF_LUT_Noise <- list()
  for (parm in Parms2Estimate){
    BRF_LUT_Noise[[parm]] <- BRF_LUT[Bands2Select[[parm]],]+BRF_LUT[Bands2Select[[parm]],]*matrix(rnorm(nrow(BRF_LUT[Bands2Select[[parm]],])*ncol(BRF_LUT[Bands2Select[[parm]],]),
                                                                                                        0,NoiseLevel[[parm]]),nrow = nrow(BRF_LUT[Bands2Select[[parm]],]))
  }

  #########################################################################
  ###                     PERFORM HYBRID INVERSION                      ###
  #########################################################################
  # train SVR for each variable and each run
  modelSVR = list()
  for (parm in Parms2Estimate){
    ColParm <- which(parm==names(InputPROSAIL))
    InputVar <- InputPROSAIL[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT_Noise[[parm]],InputVar,FigPlot = TRUE,nbEnsemble = nbModels,WithReplacement=Replacement)
  }
  return(modelSVR)
}


#' This function trains a suppot vector regression for a set of variables based on spectral data
#'
#' @param BRF_LUT numeric. LUT of bidirectional reflectances factors used for training
#' @param InputVar numeric. biophysical parameter corresponding to the reflectance
#' @param FigPlot Boolean. Set to TRUE if you want a scatterplot
#' @param nbEnsemble numeric. Number of individual subsets should be generated from BRF_LUT
#' @param WithReplacement Boolean. should subsets be generated with or without replacement?
#'
#' @return modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
#' @importFrom liquidSVM svmRegression
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom graphics par
#' @import dplyr
#' @import ggplot2
# @' @import caret
#' @export

#trainSVR_Hybrid <- function(BRF_LUT,InputVar,FigPlot = FALSE,nbEnsemble = 20,WithReplacement=FALSE){
PROSAIL_Hybrid_Train <- function(BRF_LUT,InputVar,FigPlot = FALSE,nbEnsemble = 20,WithReplacement=FALSE){

  x <- y <- ymean <- ystdmin <- ystdmax <- NULL
  # library(dplyr)
  # split the LUT into nbEnsemble subsets
  nbSamples <- length(InputVar)
  if (dim(BRF_LUT)[2]==nbSamples){
    BRF_LUT <- t(BRF_LUT)
  }

  # if subsets are generated from BRF_LUT with replacement
  if (WithReplacement==TRUE){
    Subsets <- list()
    samples_per_run <- round(nbSamples/nbEnsemble)
    for (run in (1:nbEnsemble)){
      Subsets[[run]] <- sample(seq(1,nbSamples), samples_per_run, replace = TRUE)
    }
  # if subsets are generated from BRF_LUT without replacement
  } else if (WithReplacement==FALSE){
    Subsets <- split(sample(seq(1,nbSamples,by = 1)),seq(1,nbEnsemble,by = 1))
  }

  # run training for each subset
  modelsSVR <- list()
  predictedYAll <- list()
  tunedModelYAll <- list()
  pb <- progress_bar$new(
    format = "Training SVR on subsets [:bar] :percent in :elapsed",
    total = nbEnsemble, clear = FALSE, width= 100)
  for (i in 1:nbEnsemble){
    pb$tick()
    Sys.sleep(1 / 100)
    TrainingSet <- list()
    TrainingSet$X <- BRF_LUT[Subsets[i][[1]],]
    TrainingSet$Y <- InputVar[Subsets[i][[1]]]
    # liquidSVM
    tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y)
    modelsSVR[[i]] <- tunedModel
  }

  # if scatterplots needed
  if (FigPlot==TRUE){
    # predict for full BRF_LUT
    for (i in 1:nbEnsemble){
      tunedModelY <- stats::predict(modelsSVR[[i]], BRF_LUT)
      tunedModelYAll = cbind(tunedModelYAll,matrix(tunedModelY,ncol = 1))
    }
    # plot prediction
    df <- data.frame(x = rep(1:nbSamples,nbEnsemble), y = as.numeric(matrix(tunedModelYAll,ncol = 1)))
    df.summary <- df %>% dplyr::group_by(x) %>%
      summarize( ymin = min(y),ystdmin = mean(y)-sd(y),
                 ymax = max(y),ystdmax = mean(y)+sd(y),
                 ymean = mean(y))
    par(mar=rep(.1, 4))
    p <- ggplot(df.summary, aes(x = InputVar, y = ymean)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = ystdmin, ymax = ystdmax))
    MeanPredict <- rowMeans(matrix(as.numeric(tunedModelYAll),ncol = nbEnsemble))
    print(p)
  }
  return(modelsSVR)
}

#' This function applies the regression models trained with PROSAIL_Hybrid_Train
#'
#' @param RegressionModels list. List of regression models produced by PROSAIL_Hybrid_Train
#' @param Refl numeric. LUT of bidirectional reflectances factors used for training
#'
#' @return HybridRes list. Estimated values corresponding to Refl. Includes
#' - MeanEstimate = mean value for the ensemble regression model
#' - StdEstimate = std value for the ensemble regression model
#' @importFrom stats predict
#' @importFrom matrixStats rowSds
#' @importFrom progress progress_bar
#' @export

PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl){

  # make sure Refl is right dimensions
  Refl <- t(Refl)
  nbFeatures <- RegressionModels[[1]]$dim
  if (!ncol(Refl)==nbFeatures & nrow(Refl)==nbFeatures){
    Refl <- t(Refl)
  }
  nbEnsemble <- length( RegressionModels)
  EstimatedVal <- list()
  pb <- progress_bar$new(
    format = "Applying SVR models [:bar] :percent in :elapsed",
    total = nbEnsemble, clear = FALSE, width= 100)
  for (i in 1:nbEnsemble){
    pb$tick()
    Sys.sleep(1 / 100)
    EstimatedVal[[i]] <- predict(RegressionModels[[i]], Refl)
  }
  EstimatedVal <- do.call(cbind,EstimatedVal)
  MeanEstimate <- rowMeans(EstimatedVal)
  StdEstimate <- rowSds(EstimatedVal)
  HybridRes <- list("MeanEstimate" = MeanEstimate,"StdEstimate" =StdEstimate)
  return(HybridRes)
}

#
# write_raster <- function(Image, HDR, ImagePath,parm) {
#
#   # Write image with resolution corresponding to window_size
#   HDR_output <- HDR
#   HDR_output$bands <- 1
#   HDR_output$`data type` <- 4
#   HDR_output$`band names` <- parm
#   Image_Format <- ENVI_type2bytes(HDR_output)
#   headerFpath <- paste(ImagePath, ".hdr", sep = "")
#   write_ENVI_header(HDR_output, headerFpath)
#   Image_Format <- ENVI_type2bytes(HDR_output)
#   ImgWrite <- array(Image, c(HDR_output$lines, HDR_output$samples, 1))
#   # ImgWrite <- aperm(ImgWrite, c(2, 3, 1))
#   fidOUT <- file(
#     description = ImagePath, open = "wb", blocking = TRUE,
#     encoding = getOption("encoding"), raw = FALSE
#   )
#   writeBin(c(ImgWrite), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
#   close(fidOUT)
#   return("")
# }

