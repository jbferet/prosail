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

#' This function trains a suppot vector regression for a set of variables based on spectral data
#'
#' @param BRF_LUT numeric. LUT of bidirectional reflectances factors used for training
#' @param InputVar numeric. biophysical parameter corresponding to the reflectance
#' @param FigPlot Boolean. Set to TRUE if you want a scatterplot
#' @param nbEnsemble numeric. Number of individual subsets should be generated from BRF_LUT
#' @param WithReplacement Boolean. should subsets be generated with or without replacement?
#'
#' @return modelsSVR list. regression models trained for the retroeval of InputVar based on BRF_LUT
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
      summarize(ymin = min(y),ystdmin = mean(y)-sd(y),
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

