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

#' This function applies noise defined in S2 ATBD to BRF LUT
#'
#' @param BRF_LUT numeric. BRF LUT
#'
#' @return BRF_LUT_Noise numeric.
#' @export

apply_noise_atbd <- function(BRF_LUT){
  AD <- AI <- 0.01
  MD <- MI <- 0.02
  WL_B2_B3 <- which(row.names(BRF_LUT)=='B2' | row.names(BRF_LUT)=='B3')
  WL_misc <- which(!row.names(BRF_LUT)=='B2' & !row.names(BRF_LUT)=='B3')
  BRF_LUT_Noise <- 0*BRF_LUT

  # add multiplicative noise to B2 and B3
  ADfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,AD), nrow = length(WL_B2_B3))
  AIfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,AI), nrow = length(WL_B2_B3))
  MDfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,MD), nrow = length(WL_B2_B3))
  MIfull <- matrix(rnorm(length(WL_B2_B3)*ncol(BRF_LUT),0,MI), nrow = length(WL_B2_B3))
  BRF_LUT_Noise[WL_B2_B3,] <- BRF_LUT[WL_B2_B3,]*(1+(MDfull+MIfull)) + ADfull + AIfull

  # add multiplicative noise to other bands
  ADfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,AD), nrow = length(WL_misc))
  AIfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,AI), nrow = length(WL_misc))
  MDfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,MD), nrow = length(WL_misc))
  MIfull <- matrix(rnorm(length(WL_misc)*ncol(BRF_LUT),0,MI), nrow = length(WL_misc))
  BRF_LUT_Noise[WL_misc,] <- BRF_LUT[WL_misc,]*(1+MDfull) + ADfull
  return(BRF_LUT_Noise)
}


#' This function applies SVR model on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param HybridModel list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a vegetation parameter
#' @param PathOut character. path for directory where results are written
#' @param SelectedBands list. list of spectral bands to be selected from raster (identified by name of vegetation parameter)
#' @param bandname character. spectral bands corresponding to the raster
#' @param MaskRaster character. path for binary mask defining ON (1) and OFF (0) pixels in the raster
#' @param MultiplyingFactor numeric. multiplying factor used to write reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected in the same range
#' @param method character. which machine learning regression method should be used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come
#'
#' @return res character. path for output files corresponding to biophysical properties
#' @importFrom progress progress_bar
#' @importFrom stars read_stars
#' @importFrom raster raster brick blockSize readStart readStop getValues writeStart writeStop writeValues
#' @importFrom matrixStats rowSds
#' @export
Apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
                                    SelectedBands, bandname,
                                    MaskRaster = FALSE,MultiplyingFactor=10000,
                                    method = 'liquidSVM'){

  # explain which biophysical variables will be computed
  BPvar <- names(HybridModel)
  # Produce a map for each biophysical property
  BPvarpath <- BPvarSDpath <- list()
  for (parm in BPvar){
    print(paste('Computing',parm,sep = ' '))
    # read by chunk to avoid memory problem
    blk <- blockSize(brick(raster_path))
    # reflectance file
    r_in <- readStart(brick(raster_path))
    # mask file
    r_inmask <- FALSE
    if (MaskRaster==FALSE){
      SelectPixels <- 'ALL'
    } else if (!MaskRaster==FALSE){
      if (file.exists(MaskRaster)){
        r_inmask <- readStart(raster(MaskRaster))
      } else if (!file.exists(MaskRaster)){
        message('WARNING: Mask file does not exist:')
        print(MaskRaster)
        message('Processing all image')
        SelectPixels <- 'ALL'
      }
    }
    # initiate progress bar
    pgbarlength <- length(HybridModel[[parm]])*blk$n
    pb <- progress_bar$new(
      format = "Hybrid inversion on raster [:bar] :percent in :elapsedfull , estimated time remaining :eta",
      total = pgbarlength, clear = FALSE, width= 100)

    # output files
    BPvarpath[[parm]] <- file.path(PathOut,paste(basename(raster_path),parm,sep = '_'))
    BPvarSDpath[[parm]] <- file.path(PathOut,paste(basename(raster_path),parm,'STD',sep = '_'))
    r_outMean <- writeStart(raster(raster_path), filename = BPvarpath[[parm]],format = "ENVI", overwrite = TRUE)
    r_outSD <- writeStart(raster(raster_path), filename = BPvarSDpath[[parm]],format = "ENVI", overwrite = TRUE)
    Selbands <- match(SelectedBands[[parm]],bandname)

    # loop over blocks
    for (i in seq_along(blk$row)) {
      # read values for block
      # format is a matrix with rows the cells values and columns the layers
      BlockVal <- getValues(r_in, row = blk$row[i], nrows = blk$nrows[i])
      FullLength <- dim(BlockVal)[1]

      if (typeof(r_inmask)=='logical'){
        BlockVal <- BlockVal[,Selbands]
        # automatically filter pixels corresponding to negative values
        SelectPixels <- which(BlockVal[,1]>0)
        BlockVal <- BlockVal[SelectPixels,]
      } else if (typeof(r_inmask)=='S4'){
        MaskVal <- getValues(r_inmask, row = blk$row[i], nrows = blk$nrows[i])
        SelectPixels <- which(MaskVal ==1)
        BlockVal <- BlockVal[SelectPixels,Selbands]
      }
      # add name for variables
      if (!method =='liquidSVM'){
        colnames(BlockVal) <- colnames(HybridModel[[parm]][[1]]$trainingData)[-1]
      }

      Mean_EstimateFull <- NA*vector(length = FullLength)
      STD_EstimateFull <- NA*vector(length = FullLength)
      if (length(SelectPixels)>0){
        BlockVal <- BlockVal/MultiplyingFactor
        modelSVR_Estimate <- list()
        for (modind in 1:length(HybridModel[[parm]])){
          pb$tick()
          modelSVR_Estimate[[modind]] <- predict(HybridModel[[parm]][[modind]], BlockVal)
        }
        modelSVR_Estimate <- do.call(cbind,modelSVR_Estimate)
        # final estimated value = mean parm value for all models
        Mean_Estimate <- rowMeans(modelSVR_Estimate)
        # 'uncertainty' = STD value for all models
        STD_Estimate <- matrixStats::rowSds(modelSVR_Estimate)
        Mean_EstimateFull[SelectPixels] <- Mean_Estimate
        STD_EstimateFull[SelectPixels] <- STD_Estimate
      } else {
        for (modind in 1:length(HybridModel[[parm]])){
          pb$tick()
        }
      }
      r_outMean <- writeValues(r_outMean, Mean_EstimateFull, blk$row[i],format = "ENVI", overwrite = TRUE)
      r_outSD <- writeValues(r_outSD, STD_EstimateFull, blk$row[i],format = "ENVI", overwrite = TRUE)
    }
    # close files
    r_in <- readStop(r_in)
    if (typeof(r_inmask)=='S4'){
      r_inmask <- readStop(r_inmask)
    }
    r_outMean <- writeStop(r_outMean)
    r_outSD <- writeStop(r_outSD)
    # write biophysical variable name in headers
    HDR <- read_ENVI_header(get_HDR_name(BPvarpath[[parm]]))
    HDR$`band names` <- paste('{',parm,'}',sep = '')
    write_ENVI_header(HDR, get_HDR_name(BPvarpath[[parm]]))
  }
  print('processing completed')
  res <- list('mean' = BPvarpath, 'SD' = BPvarSDpath)
  return(res)
}

#' get hdr name from image file name, assuming it is BIL format
#'
#' @param ImPath path of the image
#'
#' @return corresponding hdr
#' @importFrom tools file_ext file_path_sans_ext
#' @export
get_HDR_name <- function(ImPath) {
  if (tools::file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (tools::file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (tools::file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(tools::file_path_sans_ext(ImPath), ".hdr", sep = "")
  }

  if (!file.exists(ImPathHDR)) {
    message("WARNING : COULD NOT FIND HDR FILE")
    print(ImPathHDR)
    message("Process may stop")
  }
  return(ImPathHDR)
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
#' @importFrom methods is
#' @export

PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl){

  # make sure Refl is right dimensions
  Refl <- t(Refl)
  if (is(RegressionModels[[1]])[1]=='liquidSVM'){
    nbFeatures <- RegressionModels[[1]]$dim
  } else {
    nbFeatures <- ncol(RegressionModels[[1]]$trainingData) - 1
  }
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
    EstimatedVal[[i]] <- predict(RegressionModels[[i]], Refl)
  }
  EstimatedVal <- do.call(cbind,EstimatedVal)
  MeanEstimate <- rowMeans(EstimatedVal)
  StdEstimate <- rowSds(EstimatedVal)
  HybridRes <- list("MeanEstimate" = MeanEstimate,"StdEstimate" =StdEstimate)
  return(HybridRes)
}

#' This function trains a support vector regression for a set of variables based on spectral data
#'
#' @param BRF_LUT numeric. LUT of bidirectional reflectances factors used for training
#' @param InputVar numeric. biophysical parameter corresponding to the reflectance
#' @param FigPlot Boolean. Set to TRUE if you want a scatterplot
#' @param nbEnsemble numeric. Number of individual subsets should be generated from BRF_LUT
#' @param WithReplacement Boolean. should subsets be generated with or without replacement?
#' @param method character. which machine learning regression method should be used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training
#'
#' @return modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
#' @importFrom liquidSVM svmRegression
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom graphics par
#' @importFrom expandFunctions reset.warnings
#' @importFrom stringr str_split
#' @importFrom simsalapar tryCatch.W.E
#' @importFrom caret train trainControl
#' @import dplyr
#' @import ggplot2
#' @export

PROSAIL_Hybrid_Train <- function(BRF_LUT, InputVar, FigPlot = FALSE, nbEnsemble = 20,
                                 WithReplacement = FALSE, method = 'liquidSVM', verbose = FALSE){

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
  modelsSVR <- predictedYAll <- tunedModelYAll <- list()
  pb <- progress_bar$new(
    format = "Training SVR on subsets [:bar] :percent in :elapsedfull , eta = :eta",
    total = nbEnsemble, clear = FALSE, width= 100)
  for (i in 1:nbEnsemble){
    TrainingSet <- list()
    TrainingSet$X <- BRF_LUT[Subsets[i][[1]],]
    TrainingSet$Y <- InputVar[Subsets[i][[1]]]
    if (method == 'liquidSVM'){
      # liquidSVM
      r1 <- tryCatch.W.E(tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y))
      # reset.warnings()
      # tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y)
      if (!is.null(r1$warning)){
        Msg <- r1$warning$message
        ValGamma <- str_split(string = Msg,pattern = 'gamma=')[[1]][2]
        ValLambda <- str_split(string = Msg,pattern = 'lambda=')[[1]][2]
        if (!is.na(as.numeric(ValGamma))){
          if (verbose==T) { message('Adjusting Gamma accordingly')}
          ValGamma <- as.numeric(ValGamma)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y,min_gamma = ValGamma)
        }
        if (!is.na(as.numeric(ValLambda))){
          if (verbose==T) { message('Adjusting Lambda accordingly')}
          ValLambda <- as.numeric(ValLambda)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y,min_lambda = ValLambda)
        }
      }
    } else {
      ctrl <- caret::trainControl(
        method = "cv",
        number = 5,
      )
      if (is.null(colnames(TrainingSet$X))){
        colnames(TrainingSet$X) <- paste('var',seq(1,ncol(TrainingSet$X)),sep = '_')
      }
      target <- matrix(TrainingSet$Y,ncol = 1)
      if (is.null(colnames(target))){
        colnames(target) <- 'target'
      }
      TrainingData <- cbind(target,TrainingSet$X)
      if (method =='svmRadial'){
        tuneGrid <- expand.grid(
          C = exp(seq(-10,0,by=1)),
          sigma = exp(seq(-10,0,by=1))
        )
      } else if (method =='svmLinear'){
        tuneGrid <- expand.grid(
          C = exp(seq(-10,0,by=1))
        )
      }
      tunedModel <- caret::train(target ~ .,
                                 data = TrainingData,
                                 method = method,
                                 preProcess = c("center", "scale"),
                                 trControl = ctrl,
                                 tuneGrid = tuneGrid)
    }
    modelsSVR[[i]] <- tunedModel
    pb$tick()
  }
  # if scatterplots needed
  if (FigPlot==TRUE){
    # predict for full BRF_LUT
    if (is.null(colnames(BRF_LUT))){
      colnames(BRF_LUT) <- paste('var',seq(1,ncol(BRF_LUT)),sep = '_')
    }
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

#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export
read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath)) {
    stop("File extension should be .hdr")
  }
  HDR <- readLines(HDRpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", HDR[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    HDR <- HDR [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  HDR <- gsub("\\{([^}]*)\\}", "\\1", HDR)
  l <- grep("\\{", HDR)
  r <- grep("\\}", HDR)

  if (length(l) != length(r)) {
    stop("Error matching curly braces in header (differing numbers).")
  }

  if (any(r <= l)) {
    stop("Mismatch of curly braces in header.")
  }

  HDR[l] <- sub("\\{", "", HDR[l])
  HDR[r] <- sub("\\}", "", HDR[r])

  for (i in rev(seq_along(l))) {
    HDR <- c(
      HDR [seq_len(l [i] - 1)],
      paste(HDR [l [i]:r [i]], collapse = "\n"),
      HDR [-seq_len(r [i])]
    )
  }

  ## split key = value constructs into list with keys as names
  HDR <- sapply(HDR, split_line, "=", USE.NAMES = FALSE)
  names(HDR) <- tolower(names(HDR))

  ## process numeric values
  tmp <- names(HDR) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  HDR [tmp] <- lapply(HDR [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })

  return(HDR)
}

#' ENVI functions
#'
#' based on https://github.com/cran/hyperSpec/blob/master/R/read.ENVI.R
#' added wavelength, fwhm, ... to header reading
#' Title
#'
#' @param x character.
#' @param separator character
#' @param trim.blank boolean.
#'
#' @return list.
#' @export
split_line <- function(x, separator, trim.blank = TRUE) {
  tmp <- regexpr(separator, x)
  key <- substr(x, 1, tmp - 1)
  value <- substr(x, tmp + 1, nchar(x))
  if (trim.blank) {
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    key <- sub(blank.pattern, "\\1", key)
    value <- sub(blank.pattern, "\\1", value)
  }
  value <- as.list(value)
  names(value) <- key
  return(value)
}

#' This function performs full training for hybrid invrsion using SVR with
#' values for default parameters
#'
#' @param InputPROSAIL list. user-defined list of input parameters to be used to produce a training LUT
#' @param atbd boolean. should input parameter distribution from ATBD be applied ?
#' @param GeomAcq list. geometry of acquisiton. list should contain min and max values for tts, tto and psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#' @param minval list. minimum value for input parameters sampled to produce a training LUT
#' @param maxval list. maximum value for input parameters sampled to produce a training LUT
#' @param TypeDistrib  list. Type of distribution. Either 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean value and STD corresponding to the parameters sampled with gaussian distribution
#' @param ParmSet list. list of input parameters set to a specific value
#' @param SAILversion character. Either 4SAIL or 4SAIL2
#' @param nbSamples numeric. number of samples in training LUT
#' @param nbSamplesPerRun numeric. number of training sample per individual regression model
#' @param nbModels numeric. number of individual models to be run for ensemble
#' @param Replacement bolean. is there replacement in subsampling?
#' @param Parms2Estimate list. list of input parameters to be estimated
#' @param Bands2Select list. list of bands used for regression for each input parameter
#' @param NoiseLevel list. list of noise value added to reflectance (defined per input parm)
#' @param SRF list. Spectral response function
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for clear conditions
#' @param Path_Results character. path for results
#' @param FigPlot boolean. Set TRUE to get scatterplot of estimated biophysical variable during training step
#' @param Force4LowLAI boolean. Set TRUE to artificially reduce leaf chemical constituent content for low LAI
#' @param method character. which machine learning regression method should be used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training
#'
#' @return modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
#' @export

train_prosail_inversion <- function(InputPROSAIL = NULL,
                                    atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                                    minval = NULL, maxval = NULL,
                                    TypeDistrib = NULL, GaussianDistrib = NULL,
                                    ParmSet = NULL, SAILversion = '4SAIL',
                                    nbSamples = 2000, nbSamplesPerRun = 100, nbModels = 20,
                                    Replacement = TRUE,
                                    Parms2Estimate = 'lai', Bands2Select = NULL, NoiseLevel = NULL,
                                    SRF = NULL,
                                    SpecPROSPECT = NULL, SpecSOIL = NULL, SpecATM = NULL,
                                    Path_Results = './', FigPlot = FALSE, Force4LowLAI = TRUE,
                                    method = 'liquidSVM', verbose = FALSE){

  # default parameter values
  defaultVal <- data.frame('CHL'=40, 'CAR'=10, 'ANT' = 0, 'EWT' = 0.01, 'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5,
                           'psoil' = 0.5, 'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1,
                           'tto' = 0, 'tts' = 30, 'psi' = 80)
  # fixed parameters
  if (is.null(ParmSet)) ParmSet <- data.frame('TypeLidf' = 2, 'alpha' = 40)

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###     1- DEFINE THE LUT USED TO TRAIN THE HYBRID INVERSION          ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ##############################################################################
  #                     user-defined set of input parameters                  ##
  ##############################################################################
  if (!is.null(InputPROSAIL)){
    InputPROSAIL <- data.frame(InputPROSAIL)
    if (atbd == FALSE | !is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('parameters to generate BRF LUT provided by user in "InputPROSAIL"')
        message('following input variables will be ignored: "atbd" "minval" "maxval" "TypeDistrib" "GaussianDistrib"')
      }
    }
    # check if all parameters defined & use default value for undefined parameters
    ListParms <- c('CHL', 'CAR', 'ANT', 'EWT', 'LMA', 'BROWN', 'N', 'psoil', 'LIDFa', 'lai', 'q', 'tto', 'tts', 'psi')
    UndefinedParms <- ListParms[which(is.na(match(ListParms, names(InputPROSAIL))))]
    for (parm in UndefinedParms){
      InputPROSAIL[[parm]] <- defaultVal[[parm]]
    }
    # Set parameters
    if (length(ParmSet)>0){
      for (parm in names(ParmSet)){
        if (is.null(InputPROSAIL[[parm]])) InputPROSAIL[[parm]] <- ParmSet[[parm]]
      }
    }

  ##############################################################################
  #                         distribution defined in S2 ATBD                   ##
  ##############################################################################
  } else if (atbd==TRUE){
    if (!is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('using PROSAIL parameter distribution defined in S2 ATBD')
        message('http://step.esa.int/docs/extra/ATBD_S2ToolBox_V2.1.pdf')
        message('following input variables will be ignored: "minval" "maxval" "TypeDistrib" "GaussianDistrib"')
      }
    }
    InputPROSAIL <- get_atbd_LUT_input(nbSamples = nbSamples,
                                       GeomAcq = GeomAcq,
                                       Codist_LAI = Codist_LAI)
    # Set parameters
    if (length(ParmSet)>0){
      for (parm in names(ParmSet)){
        if (is.null(InputPROSAIL[[parm]])) InputPROSAIL[[parm]] <- ParmSet[[parm]]
      }
    }
  ##############################################################################
  #                     user defined range and distribution                   ##
  ##############################################################################
  } else {
    # define distribution for parameters to be sampled
    res <- get_default_LUT_input(TypeDistrib = TypeDistrib, GaussianDistrib = GaussianDistrib,
                                 minval = minval, maxval = maxval)
    TypeDistrib <- res$TypeDistrib
    GaussianDistrib <- res$GaussianDistrib
    minval <- res$minval
    maxval <- res$maxval
    # produce input parameters distribution
    if (SAILversion=='4SAIL'){
      InputPROSAIL <- get_distribution_input_prosail(minval, maxval, ParmSet, nbSamples,
                                                     TypeDistrib = TypeDistrib,
                                                     Mean = GaussianDistrib$Mean, Std = GaussianDistrib$Std,
                                                     Force4LowLAI = Force4LowLAI)
    } else if (SAILversion=='4SAIL2'){
      InputPROSAIL <- get_distribution_input_prosail2(minval, maxval, ParmSet, nbSamples,
                                                      TypeDistrib = TypeDistrib,
                                                      Mean = GaussianDistrib$Mean, Std = GaussianDistrib$Std,
                                                      Force4LowLAI = Force4LowLAI)
    }
    if (SAILversion=='4SAIL2'){
      # Definition of Cv & update LAI
      # vegetation covers 100% surface when > MaxLAI (Cv=1)
      MaxLAI <- min(c(maxval$lai),4)
      InputPROSAIL$Cv <- NA*InputPROSAIL$lai
      InputPROSAIL$Cv[which(InputPROSAIL$lai>MaxLAI)] <- 1
      # Cv proportional to LAI when LAI < MaxLAI
      # needs to be refined as it is very ugly, undocumented and unclear
      InputPROSAIL$Cv[which(InputPROSAIL$lai<=MaxLAI)] <- (1/MaxLAI)+InputPROSAIL$lai[which(InputPROSAIL$lai<=MaxLAI)]/(MaxLAI+1)
      InputPROSAIL$Cv <- InputPROSAIL$Cv*matrix(rnorm(length(InputPROSAIL$Cv),mean = 1,sd = 0.1))
      InputPROSAIL$Cv[which(InputPROSAIL$Cv<0)] <- 0
      InputPROSAIL$Cv[which(InputPROSAIL$Cv>1)] <- 1
      InputPROSAIL$Cv[which(InputPROSAIL$lai>MaxLAI)] <- 1
      InputPROSAIL$lai <- InputPROSAIL$lai*InputPROSAIL$Cv
      # fraction of brown vegetation is not included
      InputPROSAIL$fraction_brown <- 0+0*InputPROSAIL$lai
      # dissociation between main vegetation and brown vegetation = 0 as fraction_brown = 0
      InputPROSAIL$diss <- 0+0*InputPROSAIL$lai
      # Zeta = ratio of crown diameter to crown height
      InputPROSAIL$Zeta <- 0.2+0*InputPROSAIL$lai
    }
  }

  # define default SpecPROSPECT, SpecSOIL and SpecATM if undefined
  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prosail::SpecPROSPECT
  if (is.null(SpecSOIL)) SpecSOIL <- prosail::SpecSOIL
  if (is.null(SpecATM)) SpecATM <- prosail::SpecATM

  # generate LUT of BRF corresponding to InputPROSAIL, for a sensor
  BRF_LUT <- Generate_LUT_BRF(SAILversion = SAILversion, InputPROSAIL = InputPROSAIL,
                              SpecPROSPECT = SpecPROSPECT, SpecSOIL = SpecSOIL, SpecATM = SpecATM)

  # apply sensor spectral response function if provided
  wvl <- SpecPROSPECT$lambda
  if (!is.null(SRF)) {
    if (!length(SRF$Spectral_Bands)==nrow(BRF_LUT)){
      BRF_LUT <- applySensorCharacteristics(wvl = wvl, SRF = SRF, InRefl = BRF_LUT)
      SpecSensor <- PrepareSensorSimulation(SpecPROSPECT,SpecSOIL,SpecATM,SRF)
    }
    rownames(BRF_LUT) <- SRF$Spectral_Bands
  }

  # write parameters LUT
  output <- matrix(unlist(InputPROSAIL), ncol = length(InputPROSAIL), byrow = FALSE)
  filename <- file.path(Path_Results,'PROSAIL_LUT_InputParms.txt')
  write.table(x = format(output, digits=3),file = filename,append = F, quote = F,
              col.names = names(InputPROSAIL), row.names = F,sep = '\t')
  # Write BRF LUT corresponding to parameters LUT
  filename <- file.path(Path_Results,'PROSAIL_LUT_Reflectance.txt')
  write.table(x = format(t(BRF_LUT), digits=5),file = filename,append = F, quote = F,
              col.names = SpecSensor$SpecPROSPECT_Sensor$lambda, row.names = F,sep = '\t')

  # bands used for inversion
  for (parm in Parms2Estimate){
    if (is.null(Bands2Select[[parm]])) Bands2Select[[parm]] <- seq(1,length(SpecSensor$BandNames))
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     add noise to reflectance data                 ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # if NoiseLevel == NULL then use the same strategy than ATBD
  BRF_LUT_Noise <- list()
  if (is.null(NoiseLevel)){
    if (SRF$Sensor=='Sentinel_2' | SRF$Sensor=='Sentinel_2A' | SRF$Sensor=='Sentinel_2B'){
      BRF_LUT_NoiseAll <- apply_noise_atbd(BRF_LUT)
      for (parm in Parms2Estimate) BRF_LUT_Noise[[parm]] <- BRF_LUT_NoiseAll[Bands2Select[[parm]],]
    } else {
      for (parm in Parms2Estimate) {
        NoiseLevel[[parm]] <- 0.01
        subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
        BRF_LUT_Noise[[parm]] <- subsetRefl + subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),0,NoiseLevel[[parm]]),
                                                                nrow = nrow(subsetRefl))
      }
    }

  } else {
    # produce LUT with noise
    for (parm in Parms2Estimate){
      subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
      BRF_LUT_Noise[[parm]] <- subsetRefl + subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),0,NoiseLevel[[parm]]),
                                                              nrow = nrow(subsetRefl))
    }
  }


  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     PERFORM HYBRID INVERSION                      ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # train SVR for each variable and each run
  modelSVR = list()
  for (parm in Parms2Estimate){
    ColParm <- which(parm==names(InputPROSAIL))
    InputVar <- InputPROSAIL[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],InputVar = InputVar,
                                             FigPlot = FigPlot,nbEnsemble = nbModels,
                                             WithReplacement = Replacement, method = method, verbose = verbose)
  }
  return(modelSVR)
}

#' writes ENVI hdr file
#'
#' @param HDR content to be written
#' @param HDRpath Path of the hdr file
#'
#' @return None
#' @importFrom stringr str_count
#' @export
write_ENVI_header <- function(HDR, HDRpath) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && str_count(x, "\\w+") > 1)) {
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(HDR), h, sep = " = ")), con = HDRpath)
  return(invisible())
}
