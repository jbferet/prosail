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
#' @param MultiplyingFactor numeric. multiplying factor used to write
#' reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' @param maxRows numeric. number of rows to process from raster each time
#' @param bigRaster boolean. should R package bigRaster be used to apply prosail
#' inversion on raster data? check https://gitlab.com/jbferet/bigRaster for
#' additional support
#' @param progressBar boolean. should progressbar be displayed?
#' @param filetype character. driver for raster file
#'
#'
#' @return res character. path for output files corresponding to
#' biophysical properties
#' @importFrom progress progress_bar
#' @importFrom raster raster brick blockSize readStart readStop getValues
#' writeStart writeStop writeValues
#' @importFrom matrixStats rowSds
#' @importFrom tools file_path_sans_ext
#' @export

apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
                                    SelectedBands, bandname, MaskRaster = NULL,
                                    MultiplyingFactor = 10000, maxRows = 100,
                                    bigRaster = FALSE, progressBar = TRUE,
                                    filetype = 'GTiff'){
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
    funct <- bigRaster::wrapperBig_prosail_inversion
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
      if (filetype %in% c('GTiff', 'COG'))
        output_rasters[[listname]] <- paste0(output_rasters[[listname]],'.tiff')
      BPvarpath[[parm]] <- output_rasters[[listname]]

      listname <- paste0('SD_',parm)
      output_rasters[[listname]] <- file.path(PathOut, paste(raster_name, parm,
                                                             'STD', sep = '_'))
      if (filetype %in% c('GTiff', 'COG'))
        output_rasters[[listname]] <- paste0(output_rasters[[listname]],'.tiff')
      BPvarSDpath[[parm]] <- output_rasters[[listname]]
    }
    bandNames <- as.list(names(output_rasters))
    names(bandNames) <- names(output_rasters)
    bigRaster::apply_bigRaster(funct = funct,
                               input_rasters = input_rasters,
                               input_args = input_args,
                               output_rasters = output_rasters,
                               output_lyrs = 1,
                               filetype = filetype,
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
      if (filetype == 'EHdr'){
        formatfile <- 'ENVI'
      } else {
        formatfile <- filetype
      }
      BPvarpath[[parm]] <- file.path(PathOut, paste(raster_name, parm,
                                                    sep = '_'))
      BPvarSDpath[[parm]] <- file.path(PathOut,paste(raster_name, parm, 'STD',
                                                     sep = '_'))
      r_outMean <- writeStart(raster(raster_path),
                              filename = BPvarpath[[parm]],
                              format = formatfile, overwrite = TRUE)
      r_outSD <- writeStart(raster(raster_path),
                            filename = BPvarSDpath[[parm]],
                            format = formatfile, overwrite = TRUE)
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
                                 format = formatfile, overwrite = TRUE)
        r_outSD <- writeValues(r_outSD, STD_EstimateFull, blk$row[i],
                               format = formatfile, overwrite = TRUE)
      }
      # close files
      r_in <- readStop(r_in)
      if (typeof(r_inmask)=='S4') r_inmask <- readStop(r_inmask)
      r_outMean <- writeStop(r_outMean)
      r_outSD <- writeStop(r_outSD)
      # write biophysical variable name in headers
      if (filetype %in% c('EHdr', 'ENVI')){
        HDR <- read_ENVI_header(get_HDR_name(BPvarpath[[parm]]))
        HDR$`band names` <- paste('{',parm,'}',sep = '')
        write_ENVI_header(HDR, get_HDR_name(BPvarpath[[parm]]))
        HDR <- read_ENVI_header(get_HDR_name(BPvarSDpath[[parm]]))
        HDR$`band names` <- paste('{',parm,'}',sep = '')
        write_ENVI_header(HDR, get_HDR_name(BPvarSDpath[[parm]]))
        BPvarpath[[parm]] <- paste0(BPvarpath[[parm]],'.envi')
        BPvarSDpath[[parm]] <- paste0(BPvarSDpath[[parm]],'.envi')
      }
      if (filetype %in% c('GTiff', 'COG')){
        BPvarpath[[parm]] <- paste0(BPvarpath[[parm]],'.tiff')
        BPvarSDpath[[parm]] <- paste0(BPvarSDpath[[parm]],'.tiff')
      }
    }
    print('processing completed')
  }
  return(list('mean' = BPvarpath, 'SD' = BPvarSDpath))
}
