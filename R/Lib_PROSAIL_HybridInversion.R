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

#' This function applies additive and multiplicative noise to BRF data
#' - if AdditiveNoise and MultiplicativeNoise are defined with a unique value,
#' noise is homogeneous across all spectrum
#' - if AdditiveNoise and MultiplicativeNoise are the same length as the number
#' of spectral bands (rows in BRF_LUT), noise is specific to each spectral band
#'
#' @param BRF_LUT numeric. BRF LUT
#' @param AdditiveNoise numeric. additive noise (0 = 0%, 1 = 100%)
#' @param MultiplicativeNoise numeric. multiplicative noise (0 = 0%, 1 = 100%)
#'
#' @return BRF_LUT_Noise numeric.
#' @export

apply_noise_AddMult <- function(BRF_LUT, AdditiveNoise = 0.01,
                                MultiplicativeNoise = 0.02){
  nbWL <- nrow(BRF_LUT)
  nbSamples <- ncol(BRF_LUT)
  # add noise to BRF
  if (length(AdditiveNoise)==1){
    AddComp <- matrix(rnorm(nbWL*nbSamples,0,AdditiveNoise),
                      nrow = nbWL)
  } else if ((length(AdditiveNoise)==nbWL)){
    AddComp <- matrix(data = 0, nrow = nbWL, ncol = nbSamples)
    for (i in seq_len(nbWL)) AddComp[i,] <- matrix(data = rnorm(nbSamples,
                                                                mean = 0,
                                                                sd = AdditiveNoise[i]),
                                            ncol = nbSamples)
  }
  if (length(MultiplicativeNoise)==1){
    MultComp <- matrix(rnorm(nbWL*nbSamples,0,MultiplicativeNoise), nrow = nbWL)
  } else if ((length(MultiplicativeNoise)==nbWL)){
    MultComp <- matrix(data = 0, nrow = nbWL, ncol = nbSamples)
    for (i in seq_len(nbWL)) MultComp[i,] <- matrix(data = rnorm(nbSamples,
                                                          mean = 0,
                                                          sd = MultiplicativeNoise[i]),
                                             ncol = nbSamples)
  }
  BRF_LUT_Noise <- BRF_LUT*(1+(MultComp)) + AddComp
  return(BRF_LUT_Noise)
}


#' This function applies SVR model on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param HybridModel list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param PathOut character. path for directory where results are written
#' @param SelectedBands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param bandname character. spectral bands corresponding to the raster
#' @param MaskRaster character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param MultiplyingFactor numeric. multiplying factor used to write reflectance
#' in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' @param maxRows numeric. number of rows to process from raster each time
#' @param bigRaster boolean. should R package bigRaster be used to apply prosail
#' inversion on raster data? check https://gitlab.com/jbferet/bigRaster for
#' additional support
#' @param progressBar boolean. should progressbar be displayed?
#'
#'
#' @return res character. path for output files corresponding to biophysical properties
#' @importFrom progress progress_bar
#' @importFrom raster raster brick blockSize readStart readStop getValues writeStart writeStop writeValues
#' @importFrom matrixStats rowSds
#' @importFrom tools file_path_sans_ext
#' @export

Apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
                                    SelectedBands, bandname, MaskRaster = NULL,
                                    MultiplyingFactor = 10000, maxRows = 100,
                                    bigRaster = FALSE, progressBar = TRUE){
  # get raster name
  raster_name <- tools::file_path_sans_ext(basename(raster_path))
  # list of biophysical variables to compute
  BPvar <- names(HybridModel)
  # define path for maps produced for each biophysical variable
  BPvarpath <- BPvarSDpath <- list()
  # check if bigRaster supported
  is_bigRaster_available <- require("bigRaster", quietly = TRUE)
  if (!is_bigRaster_available & bigRaster){
    message('bigRaster cannot be used to apply prosail inversion')
    message('raster management will be performed with the package raster instead')
    message('this may require additional computational resource')
    message('check https://gitlab.com/jbferet/bigRaster for additional support')
    bigRaster <- FALSE
  }
  if (bigRaster){
    funct <- wrapperBig_prosail_inversion
    input_args <- list('HybridModel' = HybridModel,
                       'SelectedBands' = SelectedBands,
                       'bandname' = bandname,
                       'ReflFactor' = MultiplyingFactor)
    if (inherits(x = HybridModel[[1]][[1]], what = 'liquidSVM')) {
      input_args$method <- 'liquidSVM'
    } else {
      input_args$method <- 'caret'
    }
    input_rasters <- as.list(raster_path)
    names(input_rasters)[1] <- 'img'
    if (!is.null(MaskRaster)) input_rasters$mask <- MaskRaster
    output_rasters <- list()
    for (parm in BPvar){
      listname <- paste0('Mean_',parm)
      output_rasters[[listname]] <- file.path(PathOut, paste(raster_name, parm,
                                                             sep = '_'))
      BPvarpath[[parm]] <- output_rasters[[listname]]

      listname <- paste0('SD_',parm)
      output_rasters[[listname]] <- file.path(PathOut, paste(raster_name, parm,
                                                             'STD', sep = '_'))
      BPvarSDpath[[parm]] <- output_rasters[[listname]]
    }
    bandNames <- as.list(names(output_rasters))
    names(bandNames) <- names(output_rasters)
    bigRaster::apply_bigRaster(funct = funct,
                               input_rasters = input_rasters,
                               input_args = input_args,
                               output_rasters = output_rasters,
                               output_lyrs = 1,
                               filetype = 'EHdr',
                               bandNames = bandNames,
                               maxRows = maxRows)
  } else {
    for (parm in BPvar){
      print(paste('Computing',parm,sep = ' '))
      # read by chunk to avoid memory problem
      blk <- blockSize(brick(raster_path))
      if (blk$nrows[1]>maxRows){
        nrows_indiv <- maxRows
        blkud <- list()
        blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
        if (sum(blkud$nrows)<sum(blk$nrows)){
          blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
        }
        blkud$row <- NULL
        blkud$n <- length(blkud$nrows)
        blkud$row <- c(1,cumsum(blkud$nrows)+1)[1:blkud$n]
        blk$row <- blkud$row
        blk$nrows <- blkud$nrows
        blk$n <- blkud$n
      }

      # reflectance file
      r_in <- readStart(brick(raster_path))
      # mask file
      r_inmask <- FALSE
      if (is.null(MaskRaster)){
        SelectPixels <- 'ALL'
      } else if (!is.null(MaskRaster)){
        if (file.exists(MaskRaster)){
          r_inmask <- readStart(raster(MaskRaster))
        } else if (!file.exists(MaskRaster)){
          message('WARNING: Mask file does not exist:')
          print(MaskRaster)
          message('Processing all image')
          SelectPixels <- 'ALL'
        }
      }
      if (progressBar == TRUE){
        # initiate progress bar
        pgbarlength <- length(HybridModel[[parm]])*blk$n
        pb <- progress_bar$new(
          format = "Hybrid inversion on raster [:bar] :percent in :elapsedfull , estimated time remaining :eta",
          total = pgbarlength, clear = FALSE, width= 100)
      }
      # output files
      BPvarpath[[parm]] <- file.path(PathOut, paste(raster_name, parm, sep = '_'))
      BPvarSDpath[[parm]] <- file.path(PathOut,paste(raster_name, parm, 'STD', sep = '_'))
      r_outMean <- writeStart(raster(raster_path),
                              filename = BPvarpath[[parm]],
                              format = "ENVI", overwrite = TRUE)
      r_outSD <- writeStart(raster(raster_path),
                            filename = BPvarSDpath[[parm]],
                            format = "ENVI", overwrite = TRUE)
      Selbands <- match(SelectedBands[[parm]], bandname)

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

        if (!inherits(HybridModel[[parm]][[1]], what = 'liquidSVM')){
          colnames(BlockVal) <- colnames(HybridModel[[parm]][[1]]$trainingData)[-1]
        }

        Mean_EstimateFull <- NA*vector(length = FullLength)
        STD_EstimateFull <- NA*vector(length = FullLength)
        if (length(SelectPixels)>0){
          BlockVal <- BlockVal/MultiplyingFactor
          modelSVR_Estimate <- list()
          for (modind in seq_len(length(HybridModel[[parm]]))){
            if (progressBar == TRUE) pb$tick()
            modelSVR_Estimate[[modind]] <- predict(HybridModel[[parm]][[modind]],
                                                   BlockVal)
          }
          modelSVR_Estimate <- do.call(cbind,modelSVR_Estimate)
          # final estimated value = mean parm value for all models
          Mean_Estimate <- rowMeans(modelSVR_Estimate)
          # 'uncertainty' = STD value for all models
          STD_Estimate <- matrixStats::rowSds(modelSVR_Estimate)
          Mean_EstimateFull[SelectPixels] <- Mean_Estimate
          STD_EstimateFull[SelectPixels] <- STD_Estimate
        } else {
          for (modind in seq_len(length(HybridModel[[parm]]))) pb$tick()
        }
        r_outMean <- writeValues(r_outMean, Mean_EstimateFull, blk$row[i],
                                 format = "ENVI", overwrite = TRUE)
        r_outSD <- writeValues(r_outSD, STD_EstimateFull, blk$row[i],
                               format = "ENVI", overwrite = TRUE)
      }
      # close files
      r_in <- readStop(r_in)
      if (typeof(r_inmask)=='S4') r_inmask <- readStop(r_inmask)
      r_outMean <- writeStop(r_outMean)
      r_outSD <- writeStop(r_outSD)
      # write biophysical variable name in headers
      HDR <- read_ENVI_header(get_HDR_name(BPvarpath[[parm]]))
      HDR$`band names` <- paste('{',parm,'}',sep = '')
      write_ENVI_header(HDR, get_HDR_name(BPvarpath[[parm]]))
    }
    print('processing completed')
  }
  return(list('mean' = BPvarpath, 'SD' = BPvarSDpath))
}

# Apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
#                                     SelectedBands, bandname, MaskRaster = FALSE,
#                                     MultiplyingFactor = 10000, bigRaster = FALSE){
#
#   # list of biophysical variables to compute
#   BPvar <- names(HybridModel)
#   # define path for maps produced for each biophysical variable
#   BPvarpath <- BPvarSDpath <- list()
#   # check if bigRaster supported
#   is_bigRaster_available <- require("bigRaster", quietly = T)
#   if (!is_bigRaster_available & bigRaster){
#     message('bigRaster cannot be used to apply prosail inversion')
#     message('raster management will be performed with the package raster instead')
#     message('this may require additional computational resource')
#     message('check https://gitlab.com/jbferet/bigRaster for additional support')
#     bigRaster <- FALSE
#   }
#
#
#   if (bigRaster){
#     funct <- wrapperBig_prosail_inversion
#     input_args <- list('HybridModel' = HybridModel,
#                        'SelectedBands' = SelectedBands,
#                        'bandname' = bandname,
#                        'ReflFactor' = MultiplyingFactor)
#     if (class(HybridModel[[1]][[1]]) == 'liquidSVM') {
#       input_args$method =='liquidSVM'
#     } else {
#       input_args$method =='caret'
#     }
#     input_rasters <- as.list(raster_path)
#     names(input_rasters) <- names_S2
#     if (!is.null(MTD_LaSRC)) input_args$LaSRC <- TRUE
#     output_rasters <- list('Stack' = Refl_path)
#     bandNames <- list('Stack' = names_S2)
#
#     bigRaster::apply_bigRaster(funct = funct,
#                                input_rasters = input_rasters,
#                                input_args = input_args,
#                                output_rasters,
#                                output_lyrs = length(input_rasters),
#                                filetype = 'EHdr',
#                                bandNames = bandNames)
#
#
#
#   } else {
#     for (parm in BPvar){
#       print(paste('Computing',parm,sep = ' '))
#       # read by chunk to avoid memory problem
#       blk <- blockSize(brick(raster_path))
#       # reflectance file
#       r_in <- readStart(brick(raster_path))
#       # mask file
#       r_inmask <- FALSE
#       if (MaskRaster==FALSE){
#         SelectPixels <- 'ALL'
#       } else if (!MaskRaster==FALSE){
#         if (file.exists(MaskRaster)){
#           r_inmask <- readStart(raster(MaskRaster))
#         } else if (!file.exists(MaskRaster)){
#           message('WARNING: Mask file does not exist:')
#           print(MaskRaster)
#           message('Processing all image')
#           SelectPixels <- 'ALL'
#         }
#       }
#       # initiate progress bar
#       pgbarlength <- length(HybridModel[[parm]])*blk$n
#       pb <- progress_bar$new(
#         format = "Hybrid inversion on raster [:bar] :percent in :elapsedfull , estimated time remaining :eta",
#         total = pgbarlength, clear = FALSE, width= 100)
#
#       # output files
#       raster_name <- tools::file_path_sans_ext(basename(raster_path))
#       BPvarpath[[parm]] <- file.path(PathOut, paste(raster_name, parm, sep = '_'))
#       BPvarSDpath[[parm]] <- file.path(PathOut,paste(raster_name, parm, 'STD', sep = '_'))
#       r_outMean <- writeStart(raster(raster_path),
#                               filename = BPvarpath[[parm]],
#                               format = "ENVI", overwrite = TRUE)
#       r_outSD <- writeStart(raster(raster_path),
#                             filename = BPvarSDpath[[parm]],
#                             format = "ENVI", overwrite = TRUE)
#       Selbands <- match(SelectedBands[[parm]], bandname)
#
#       # loop over blocks
#       for (i in seq_along(blk$row)) {
#         # read values for block
#         # format is a matrix with rows the cells values and columns the layers
#         BlockVal <- getValues(r_in, row = blk$row[i], nrows = blk$nrows[i])
#         FullLength <- dim(BlockVal)[1]
#
#         if (typeof(r_inmask)=='logical'){
#           BlockVal <- BlockVal[,Selbands]
#           # automatically filter pixels corresponding to negative values
#           SelectPixels <- which(BlockVal[,1]>0)
#           BlockVal <- BlockVal[SelectPixels,]
#         } else if (typeof(r_inmask)=='S4'){
#           MaskVal <- getValues(r_inmask, row = blk$row[i], nrows = blk$nrows[i])
#           SelectPixels <- which(MaskVal ==1)
#           BlockVal <- BlockVal[SelectPixels,Selbands]
#         }
#         # add name for variables
#
#         if (!inherits(HybridModel[[parm]][[1]], what = 'liquidSVM')){
#           colnames(BlockVal) <- colnames(HybridModel[[parm]][[1]]$trainingData)[-1]
#         }
#
#         Mean_EstimateFull <- NA*vector(length = FullLength)
#         STD_EstimateFull <- NA*vector(length = FullLength)
#         if (length(SelectPixels)>0){
#           BlockVal <- BlockVal/MultiplyingFactor
#           modelSVR_Estimate <- list()
#           for (modind in seq_len(length(HybridModel[[parm]]))){
#             pb$tick()
#             modelSVR_Estimate[[modind]] <- predict(HybridModel[[parm]][[modind]],
#                                                    BlockVal)
#           }
#           modelSVR_Estimate <- do.call(cbind,modelSVR_Estimate)
#           # final estimated value = mean parm value for all models
#           Mean_Estimate <- rowMeans(modelSVR_Estimate)
#           # 'uncertainty' = STD value for all models
#           STD_Estimate <- matrixStats::rowSds(modelSVR_Estimate)
#           Mean_EstimateFull[SelectPixels] <- Mean_Estimate
#           STD_EstimateFull[SelectPixels] <- STD_Estimate
#         } else {
#           for (modind in seq_len(length(HybridModel[[parm]]))) pb$tick()
#         }
#         r_outMean <- writeValues(r_outMean, Mean_EstimateFull, blk$row[i],
#                                  format = "ENVI", overwrite = TRUE)
#         r_outSD <- writeValues(r_outSD, STD_EstimateFull, blk$row[i],
#                                format = "ENVI", overwrite = TRUE)
#       }
#       # close files
#       r_in <- readStop(r_in)
#       if (typeof(r_inmask)=='S4') r_inmask <- readStop(r_inmask)
#       r_outMean <- writeStop(r_outMean)
#       r_outSD <- writeStop(r_outSD)
#       # write biophysical variable name in headers
#       HDR <- read_ENVI_header(get_HDR_name(BPvarpath[[parm]]))
#       HDR$`band names` <- paste('{',parm,'}',sep = '')
#       write_ENVI_header(HDR, get_HDR_name(BPvarpath[[parm]]))
#     }
#     print('processing completed')
#   }
#   return(list('mean' = BPvarpath, 'SD' = BPvarSDpath))
# }

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

#' generate InputPROSAIL, following
#' - either distribution defined in ATBD
#' - or distribution defined in user
#'
#' @param atbd boolean. should input parameter distribution from ATBD be applied ?
#' @param GeomAcq list. geometry of acquisiton. list should contain min and max values for tts, tto and psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#' @param minval list. minimum value for input parameters sampled to produce a training LUT
#' @param maxval list. maximum value for input parameters sampled to produce a training LUT
#' @param TypeDistrib  list. Type of distribution. Either 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean value and STD corresponding to the parameters sampled with gaussian distribution
#' @param ParmSet list. list of input parameters set to a specific value
#' @param nbSamples numeric. number of samples in training LUT
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training
#'
#' @return InputPROSAIL
#' @export

get_InputPROSAIL <- function(atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                             minval = NULL, maxval = NULL,
                             TypeDistrib = NULL, GaussianDistrib = NULL,
                             ParmSet = NULL, nbSamples = 2000, verbose = FALSE){

  # default parameter values
  ListParms <- c('CHL', 'CAR', 'ANT', 'EWT', 'LMA', 'BROWN', 'N', 'psoil',
                 'LIDFa', 'lai', 'q', 'tto', 'tts', 'psi')
  defaultVal <- data.frame('CHL'=40, 'CAR'=10, 'ANT' = 0, 'EWT' = 0.01,
                           'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5, 'psoil' = 0.5,
                           'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1,
                           'tto' = 0, 'tts' = 30, 'psi' = 80)

  if (atbd==TRUE){
    #__________________________________________________________________________#
    #                         distribution defined in S2 ATBD                 ##
    #__________________________________________________________________________#
    if (!is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('using PROSAIL parameter distribution defined in S2 ATBD')
        message('http://step.esa.int/docs/extra/ATBD_S2ToolBox_V2.1.pdf')
      }
    }
    InputPROSAIL <- get_atbd_LUT_input(nbSamples = nbSamples,
                                       GeomAcq = GeomAcq,
                                       Codist_LAI = Codist_LAI)
  } else {
    #__________________________________________________________________________#
    #                     user defined range and distribution                 ##
    #__________________________________________________________________________#
    # check consistency between user-defined variables, use default distribution if needed
    res <- get_default_LUT_input(TypeDistrib = TypeDistrib,
                                 GaussianDistrib = GaussianDistrib,
                                 minval = minval, maxval = maxval)
    TypeDistrib <- res$TypeDistrib
    GaussianDistrib <- res$GaussianDistrib
    minval <- res$minval
    maxval <- res$maxval

    # fixed parameters
    if (is.na(match('TypeLidf', names(ParmSet)))){
      if (is.null(ParmSet)){
        ParmSet <- data.frame('TypeLidf' = 2)
      } else {
        ParmSet <- data.frame(ParmSet, 'TypeLidf' = 2)
      }
    }
    if (is.na(match('alpha', names(ParmSet)))){
      if (is.null(ParmSet)){
        ParmSet <- data.frame('alpha' = 40)
      } else {
        ParmSet <- data.frame(ParmSet, 'alpha' = 40)
      }
    }
    # produce input parameters distribution
    InputPROSAIL <- get_distribution_input_prosail(minval, maxval,
                                                   ParmSet, nbSamples,
                                                   TypeDistrib = TypeDistrib,
                                                   Mean = GaussianDistrib$Mean,
                                                   Std = GaussianDistrib$Std)
  }
  InputPROSAIL <- data.frame(InputPROSAIL)
  return(InputPROSAIL)
}

#' This function applies the regression models trained with PROSAIL_Hybrid_Train
#'
#' @param RegressionModels list. List of regression models produced by PROSAIL_Hybrid_Train
#' @param Refl numeric. LUT of bidirectional reflectances factors used for training
#' @param progressBar boolean. should progressbar be displayed?
#'
#' @return HybridRes list. Estimated values corresponding to Refl. Includes
#' - MeanEstimate = mean value for the ensemble regression model
#' - StdEstimate = std value for the ensemble regression model
#' @importFrom stats predict
#' @importFrom matrixStats rowSds
#' @importFrom progress progress_bar
#' @importFrom methods is
#' @export

PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl, progressBar = FALSE){

  # make sure Refl is right dimensions
  Refl <- t(Refl)
  if (inherits(RegressionModels[[1]], what = 'liquidSVM')){
    nbFeatures <- RegressionModels[[1]]$dim
  } else {
    nbFeatures <- ncol(RegressionModels[[1]]$trainingData) - 1
  }
  if (!ncol(Refl)==nbFeatures & nrow(Refl)==nbFeatures){
    Refl <- t(Refl)
  }
  nbEnsemble <- length( RegressionModels)
  EstimatedVal <- list()
  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Applying SVR models [:bar] :percent in :elapsed",
      total = nbEnsemble, clear = FALSE, width= 100)
  }
  for (i in seq_len(length(nbEnsemble))){
    EstimatedVal[[i]] <- predict(RegressionModels[[i]], Refl)
    if (progressBar == TRUE) pb$tick()
  }
  EstimatedVal <- do.call(cbind,EstimatedVal)
  MeanEstimate <- rowMeans(EstimatedVal)
  StdEstimate <- rowSds(EstimatedVal)
  return(list("MeanEstimate" = MeanEstimate,
              "StdEstimate" = StdEstimate))
}

#' This function trains a support vector regression for a set of variables based on spectral data
#'
#' @param BRF_LUT numeric. LUT of bidirectional reflectances factors used for training
#' @param InputVar numeric. biophysical parameter corresponding to the reflectance
#' @param nbEnsemble numeric. Number of individual subsets should be generated from BRF_LUT
#' @param WithReplacement Boolean. should subsets be generated with or without replacement?
#' @param method character. which machine learning regression method should be used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training
#' @param progressBar boolean. should progressbar be displayed?
#'
#' @return modelsMLR list. ML regression models trained for the retrieval of
#' InputVar based on BRF_LUT
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom simsalapar tryCatch.W.E
#' @importFrom caret train trainControl
#' @importFrom magrittr %>%
#' @export

PROSAIL_Hybrid_Train <- function(BRF_LUT, InputVar, nbEnsemble = 20,
                                 WithReplacement = FALSE,
                                 method = 'liquidSVM',
                                 verbose = FALSE, progressBar = FALSE){

  x <- y <- ymean <- ystdmin <- ystdmax <- NULL
  # split the LUT into nbEnsemble subsets
  nbSamples <- length(InputVar)
  if (dim(BRF_LUT)[2]==nbSamples) BRF_LUT <- t(BRF_LUT)
  # if subsets are generated from BRF_LUT with replacement
  if (WithReplacement==TRUE){
    Subsets <- list()
    samples_per_run <- round(nbSamples/nbEnsemble)
    for (run in seq_len(nbEnsemble)){
      Subsets[[run]] <- sample(seq_len(nbSamples), samples_per_run, replace = TRUE)
    }
  # if subsets are generated from BRF_LUT without replacement
  } else if (WithReplacement==FALSE){
    Subsets <- split(sample(seq_len(nbSamples)),seq_len(nbEnsemble))
  }

  # run training for each subset
  modelsMLR <- predictedYAll <- tunedModelYAll <- list()

  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Training SVR on subsets [:bar] :percent in :elapsedfull , eta = :eta",
      total = nbEnsemble, clear = FALSE, width= 100)
  }
  for (i in seq_len(nbEnsemble)){
    TrainingSet <- list()
    TrainingSet$X <- BRF_LUT %>% dplyr::slice(Subsets[i][[1]])
    TrainingSet$Y <- InputVar[Subsets[i][[1]]]

    # if using caret
    control <- caret::trainControl(method="cv", number = 5)
    if (is.null(colnames(TrainingSet$X))) colnames(TrainingSet$X) <- paste('var',seq_len(ncol(TrainingSet$X)),sep = '_')
    target <- matrix(TrainingSet$Y,ncol = 1)
    if (is.null(colnames(target))) colnames(target) <- 'target'
    TrainingData <- cbind(target,TrainingSet$X)

    if (method == 'liquidSVM'){
      # liquidSVM
      r1 <- tryCatch.W.E(tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                                TrainingSet$Y))
      # tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y)
      if (!is.null(r1$warning)){
        Msg <- r1$warning$message
        ValGamma <- stringr::str_split(string = Msg,pattern = 'gamma=')[[1]][2]
        ValLambda <- stringr::str_split(string = Msg,pattern = 'lambda=')[[1]][2]
        if (!is.na(as.numeric(ValGamma))){
          if (verbose==T) { message('Adjusting Gamma accordingly')}
          ValGamma <- as.numeric(ValGamma)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                 TrainingSet$Y,
                                                 min_gamma = ValGamma)
        }
        if (!is.na(as.numeric(ValLambda))){
          if (verbose==T) { message('Adjusting Lambda accordingly')}
          ValLambda <- as.numeric(ValLambda)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                 TrainingSet$Y,
                                                 min_lambda = ValLambda)
        }
      }
    } else if (method %in% c('svmRadial', 'svmLinear')){
      if (method =='svmRadial') tuneGrid <- expand.grid(C = exp(seq(-10, 0)),
                                                        sigma = exp(seq(-10, 0)))
      if (method =='svmLinear')tuneGrid <- expand.grid(C = exp(seq(-10, 0)))
      tunedModel <- caret::train(target ~ .,
                                 data = TrainingData,
                                 method = method,
                                 preProcess = c("center", "scale"),
                                 trControl = control,
                                 tuneGrid = tuneGrid)
    } else if (method == 'gaussprLinear'){
      tunedModel <- caret::train(target ~ .,
                                 data = TrainingData,
                                 method = method,
                                 preProcess = c("center", "scale"),
                                 trControl = control)
    } else if (method == 'rf'){
      seed <- 7
      set.seed(seed)
      mtry <- ncol(TrainingSet$X)/3
      tunedModel <- train(target~.,
                          data = TrainingData,
                          method = method,
                          metric = 'RMSE',
                          preProcess = c("center", "scale"),
                          tuneLength = 15,
                          trControl = control)
    # } else if (method == 'xgbLinear'){
    #   grid_default <- expand.grid(nrounds = 100,
    #                               max_depth = 6,
    #                               eta = 0.3,
    #                               gamma = 0,
    #                               colsample_bytree = 1,
    #                               min_child_weight = 1,
    #                               subsample = 1)
    #   train_control <- caret::trainControl(method = "none",
    #                                        verboseIter = FALSE,
    #                                        allowParallel = TRUE )
    #
    #   tunedModel <- caret::train(
    #     x = trainMNX,
    #     y = trainMNY,
    #     trControl = train_control,
    #     tuneGrid = grid_default,
    #     method = "xgbTree",
    #     verbose = TRUE,
    #     nthreads = 4
    #   )
    #
    #
    #   grid <- expand.grid(nrounds = c(10,20), lambda = c(0.1),
    #                       alpha = c(1), eta = c(0.1))
    #   tunedModel <- train(target~.,
    #                       data = TrainingData,
    #                       method = method,
    #                       metric = 'RMSE',
    #                       preProcess = c("center", "scale"),
    #                       tuneGrid = grid,
    #                       gamma = 0.5,
    #                       trControl = control)
    }
    modelsMLR[[i]] <- tunedModel
    if (progressBar == TRUE) pb$tick()
  }
  return(modelsMLR)
}

#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export

read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath) & !grepl(".HDR$", HDRpath)) {
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

  if (any(r <= l)) stop("Mismatch of curly braces in header.")

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

#' This function performs full training for hybrid inversion using SVR with
#' values for default parameters
#'
#' @param InputPROSAIL list. user-defined list of input parameters to be used to produce a training LUT
#' @param BRF_LUT list. user-defined BRF LUT used to run the hybrid inversion
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
#' @param method character. which machine learning regression method should be used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training
#'
#' @return modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
#' @export

train_prosail_inversion <- function(InputPROSAIL = NULL, BRF_LUT = NULL,
                                    atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                                    minval = NULL, maxval = NULL,
                                    TypeDistrib = NULL, GaussianDistrib = NULL,
                                    ParmSet = NULL, SAILversion = '4SAIL',
                                    nbSamples = 2000, nbSamplesPerRun = 100,
                                    nbModels = 20, Replacement = TRUE,
                                    Parms2Estimate = 'lai', Bands2Select = NULL,
                                    NoiseLevel = NULL, SRF = NULL,
                                    SpecPROSPECT = NULL, SpecSOIL = NULL,
                                    SpecATM = NULL,
                                    Path_Results = './', FigPlot = FALSE,
                                    method = 'liquidSVM', verbose = FALSE){

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###     1- DEFINE THE LUT USED TO TRAIN THE HYBRID INVERSION          ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # default parameter values
  defaultVal <- data.frame('CHL'=40, 'CAR'=10, 'ANT' = 0, 'EWT' = 0.01,
                           'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5, 'psoil' = 0.5,
                           'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1, 'tto' = 0,
                           'tts' = 30, 'psi' = 80, 'TypeLidf' = 2, 'alpha' = 40)
  ListParms <- names(defaultVal)

  ##############################################################################
  #                     user-defined set of input parameters                  ##
  ##############################################################################
  if (!is.null(InputPROSAIL)){
    InputPROSAIL <- data.frame(InputPROSAIL)
    if (atbd == TRUE | !is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('parameters to generate BRF LUT provided by user in "InputPROSAIL"')
        message('following input variables will be ignored: "atbd" "minval" "maxval" "TypeDistrib" "GaussianDistrib"')
      }
    }
    # check if all parameters defined & use default value for undefined parameters
    UndefinedParms <- ListParms[which(is.na(match(ListParms, names(InputPROSAIL))))]
    for (parm in UndefinedParms) InputPROSAIL[[parm]] <- defaultVal[[parm]]
    # Set parameters
    if (length(ParmSet)>0){
      for (parm in names(ParmSet)){
        if (is.null(InputPROSAIL[[parm]])) InputPROSAIL[[parm]] <- ParmSet[[parm]]
      }
    }
  } else {
    InputPROSAIL <- get_InputPROSAIL(atbd = atbd, GeomAcq = GeomAcq,
                                     Codist_LAI = Codist_LAI,
                                     minval = minval, maxval = maxval,
                                     TypeDistrib = TypeDistrib,
                                     GaussianDistrib = GaussianDistrib,
                                     ParmSet = ParmSet, nbSamples = nbSamples,
                                     verbose = verbose)
  }

  if (!is.null(BRF_LUT)){
    for (parm in Parms2Estimate){
      if (!is.null(BRF_LUT[[parm]])){
        BRF_LUT_Noise[[parm]] <- BRF_LUT[[parm]]
      } else {
        message('Please make sure you provide BRF_LUT as a list with elements corresponding to Parms2Estimate')
      }
    }
  } else {
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    ### 2- PRODUCE BRF from InputPROSAIL & ddefault spectral sampling = 1nm  ###
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    # define default SpecPROSPECT, SpecSOIL and SpecATM if undefined
    if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
    if (is.null(SpecSOIL)) SpecSOIL <- prosail::SpecSOIL
    if (is.null(SpecATM)) SpecATM <- prosail::SpecATM
    # check if same spectral sampling for all key variables
    check_SpectralSampling(SpecPROSPECT, SpecSOIL, SpecATM)
    # generate LUT of BRF corresponding to InputPROSAIL, for a sensor

    if (!'fCover' %in% Parms2Estimate &
        !'albedo' %in% Parms2Estimate &
        !'fAPAR' %in% Parms2Estimate){
      BRF_LUT <- Generate_LUT_BRF(SAILversion = SAILversion,
                                  InputPROSAIL = InputPROSAIL,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM)
    } else if ('fCover' %in% Parms2Estimate |
               'albedo' %in% Parms2Estimate |
               'fAPAR' %in% Parms2Estimate){
      res <- Generate_LUT_PROSAIL(SAILversion = SAILversion,
                                  InputPROSAIL = InputPROSAIL,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM)
      BRF_LUT <- res$BRF
      InputPROSAIL$fCover <- res$fCover
      InputPROSAIL$fAPAR <- res$fAPAR
      InputPROSAIL$albedo <- res$albedo
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     3- APPLY SPECTRAL RESPONSE FUNCTION if not already applied    ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # apply sensor spectral response function if provided
    wvl <- SpecPROSPECT$lambda
    if (!is.null(SRF)) {
      if (!length(SRF$Spectral_Bands)==nrow(BRF_LUT)){
        BRF_LUT <- applySensorCharacteristics(wvl = wvl, SRF = SRF, InRefl = BRF_LUT)
        SpecSensor <- PrepareSensorSimulation(SpecPROSPECT, SpecSOIL, SpecATM, SRF)
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
      if (is.null(Bands2Select[[parm]])) Bands2Select[[parm]] <- seq_len(length(SpecSensor$BandNames))
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     4- add noise to reflectance data                              ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # if NoiseLevel == NULL then use the same strategy than ATBD
    BRF_LUT_Noise <- list()
    if (is.null(NoiseLevel)){
      if (SRF$Sensor %in% c('Sentinel_2', 'Sentinel_2A', 'Sentinel_2B')){
        BRF_LUT_NoiseAll <- apply_noise_atbd(BRF_LUT)
        for (parm in Parms2Estimate) BRF_LUT_Noise[[parm]] <- BRF_LUT_NoiseAll[Bands2Select[[parm]],]
      } else {
        for (parm in Parms2Estimate) {
          NoiseLevel[[parm]] <- 0.01
          subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
          BRF_LUT_Noise[[parm]] <- subsetRefl + subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                                                        mean = 0,
                                                                        sd = NoiseLevel[[parm]]),
                                                                  nrow = nrow(subsetRefl))
        }
      }
    } else {
      # produce LUT with noise
      for (parm in Parms2Estimate){
        if (is.null(NoiseLevel[[parm]])) NoiseLevel[[parm]] <- 0.01
        subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
        BRF_LUT_Noise[[parm]] <- subsetRefl + subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                                                      mean = 0,
                                                                      sd = NoiseLevel[[parm]]),
                                                                nrow = nrow(subsetRefl))
      }
    }
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     PERFORM HYBRID INVERSION                      ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # train SVR for each variable and each run
  modelSVR <- list()
  for (parm in Parms2Estimate){
    ColParm <- which(parm==names(InputPROSAIL))
    InputVar <- InputPROSAIL[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                             InputVar = InputVar,
                                             nbEnsemble = nbModels,
                                             WithReplacement = Replacement,
                                             method = method, verbose = verbose)
  }
  return(modelSVR)
}

#' writes ENVI hdr file
#'
#' @param HDR content to be written
#' @param HDRpath Path of the hdr file
#'
#' @return None
#' @export

write_ENVI_header <- function(HDR, HDRpath) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && stringr::str_count(x, "\\w+") > 1)) {
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(HDR), h, sep = " = ")), con = HDRpath)
  return(invisible())
}
