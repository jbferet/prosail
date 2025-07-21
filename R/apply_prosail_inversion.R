#' This function applies SVR model on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param hybrid_model list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param output_path character. path for directory where results are written
#' @param selected_bands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param bandname character. spectral bands corresponding to the raster
#' @param mask_path character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param multiplying_factor numeric. multiplying factor used to write
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
#' @importFrom utils installed.packages
#' @export

apply_prosail_inversion <- function(raster_path, hybrid_model, output_path,
                                    selected_bands, bandname, mask_path = NULL,
                                    multiplying_factor = 10000, maxRows = 100,
                                    bigRaster = FALSE, progressBar = TRUE,
                                    filetype = 'GTiff'){

  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  # get raster name
  raster_name <- tools::file_path_sans_ext(basename(raster_path))
  # list of biophysical variables to compute
  bp_vars <- names(hybrid_model)
  # define path for maps produced for each biophysical variable
  bp_var_path <- bp_var_sd_path <- list()
  # check if bigRaster supported
  is_bigRaster_available <- ("bigRaster" %in% rownames(installed.packages()))
  if (!is_bigRaster_available & bigRaster){
    message('bigRaster cannot be used to apply prosail inversion')
    message('rasters will be processed with the package terra instead')
    message('this may require additional computational resource')
    message('check https://gitlab.com/jbferet/bigRaster for additional support')
    bigRaster <- FALSE
  }
  if (bigRaster){
    funct <- bigRaster::wrapperBig_prosail_inversion
    input_args <- list('hybrid_model' = hybrid_model,
                       'selected_bands' = selected_bands,
                       'bandname' = bandname,
                       'ReflFactor' = multiplying_factor)
    if (inherits(x = hybrid_model[[1]][[1]], what = 'liquidSVM')) {
      input_args$method <- 'liquidSVM'
    } else {
      input_args$method <- 'caret'
    }
    input_rasters <- as.list(raster_path)
    names(input_rasters)[1] <- 'img'
    if (!is.null(mask_path)) input_rasters$mask <- mask_path
    output_rasters <- list()
    for (parm in bp_vars){
      listname <- paste0('Mean_',parm)
      output_rasters[[listname]] <- file.path(output_path, paste(raster_name,
                                                                 parm,
                                                                 sep = '_'))
      if (filetype %in% c('GTiff', 'COG'))
        output_rasters[[listname]] <- paste0(output_rasters[[listname]],'.tiff')
      bp_var_path[[parm]] <- output_rasters[[listname]]

      listname <- paste0('SD_',parm)
      output_rasters[[listname]] <- file.path(output_path, paste(raster_name,
                                                                 parm,
                                                                 'STD',
                                                                 sep = '_'))
      if (filetype %in% c('GTiff', 'COG'))
        output_rasters[[listname]] <- paste0(output_rasters[[listname]],'.tiff')
      bp_var_sd_path[[parm]] <- output_rasters[[listname]]
    }
    band_names <- as.list(names(output_rasters))
    names(band_names) <- names(output_rasters)
    bigRaster::apply_bigRaster(funct = funct,
                               input_rasters = input_rasters,
                               input_args = input_args,
                               output_rasters = output_rasters,
                               output_lyrs = 1,
                               filetype = filetype,
                               band_names = band_names,
                               maxRows = maxRows)
  } else {
    for (parm in bp_vars){
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
      if (is.null(mask_path)){
        select_pixels <- 'ALL'
      } else if (!is.null(mask_path)){
        if (file.exists(mask_path)){
          r_inmask <- readStart(raster(mask_path))
        } else if (!file.exists(mask_path)){
          message('WARNING: Mask file does not exist:')
          print(mask_path)
          message('Processing all image')
          select_pixels <- 'ALL'
        }
      }
      if (progressBar == TRUE){
        # initiate progress bar
        pgbarlength <- length(hybrid_model[[parm]])*blk$n
        pb <- progress_bar$new(
          format = "Hybrid inversion on raster [:bar] :percent in :elapsedfull ,
          estimated time remaining :eta",
          total = pgbarlength, clear = FALSE, width= 100)
      }
      # output files
      if (filetype == 'EHdr'){
        formatfile <- 'ENVI'
      } else {
        formatfile <- filetype
      }
      bp_var_path[[parm]] <- file.path(output_path, paste(raster_name, parm,
                                                          sep = '_'))
      bp_var_sd_path[[parm]] <- file.path(output_path,paste(raster_name, parm,
                                                            'STD', sep = '_'))
      r_out_mean <- writeStart(raster(raster_path),
                               filename = bp_var_path[[parm]],
                               format = formatfile, overwrite = TRUE)
      r_out_sd <- writeStart(raster(raster_path),
                             filename = bp_var_sd_path[[parm]],
                             format = formatfile, overwrite = TRUE)
      sel_bands <- match(selected_bands[[parm]], bandname)

      # loop over blocks
      for (i in seq_along(blk$row)) {
        # read values for block
        # format is a matrix with rows the cells values and columns the layers
        block_val <- getValues(r_in, row = blk$row[i], nrows = blk$nrows[i])
        full_length <- dim(block_val)[1]

        if (typeof(r_inmask)=='logical'){
          block_val <- block_val[,sel_bands]
          # automatically filter pixels corresponding to negative values
          select_pixels <- which(block_val[,1]>0)
          block_val <- block_val[select_pixels,]
        } else if (typeof(r_inmask)=='S4'){
          MaskVal <- getValues(r_inmask, row = blk$row[i], nrows = blk$nrows[i])
          select_pixels <- which(MaskVal ==1)
          block_val <- block_val[select_pixels,sel_bands]
        }
        # add name for variables

        if (!inherits(hybrid_model[[parm]][[1]], what = 'liquidSVM')){
          colnames(block_val) <- colnames(hybrid_model[[parm]][[1]]$trainingData)[-1]
        }

        mean_estimate_full <- NA*vector(length = full_length)
        sd_estimate_full <- NA*vector(length = full_length)
        if (length(select_pixels)>0){
          block_val <- block_val/multiplying_factor
          estimates <- list()
          for (modind in seq_len(length(hybrid_model[[parm]]))){
            if (progressBar == TRUE)
              pb$tick()
            estimates[[modind]] <- predict(hybrid_model[[parm]][[modind]],
                                           block_val)
          }
          estimates <- do.call(cbind,estimates)
          # final estimated value = mean parm value for all models
          mean_estimate <- rowMeans(estimates)
          # 'uncertainty' = STD value for all models
          sd_estimate <- matrixStats::rowSds(estimates)
          mean_estimate_full[select_pixels] <- mean_estimate
          sd_estimate_full[select_pixels] <- sd_estimate
        } else {
          for (modind in seq_len(length(hybrid_model[[parm]])))
            if (progressBar == TRUE)
              pb$tick()
        }
        r_out_mean <- writeValues(r_out_mean, mean_estimate_full, blk$row[i],
                                  format = formatfile, overwrite = TRUE)
        r_out_sd <- writeValues(r_out_sd, sd_estimate_full, blk$row[i],
                                format = formatfile, overwrite = TRUE)
      }
      # close files
      r_in <- readStop(r_in)
      if (typeof(r_inmask)=='S4')
        r_inmask <- readStop(r_inmask)
      r_out_mean <- writeStop(r_out_mean)
      r_out_sd <- writeStop(r_out_sd)
      # write biophysical variable name in headers
      if (filetype %in% c('EHdr', 'ENVI')){
        hdr <- read_envi_header(get_hdr_name(bp_var_path[[parm]]))
        hdr$`band names` <- paste('{',parm,'}',sep = '')
        write_envi_header(hdr, get_hdr_name(bp_var_path[[parm]]))
        hdr <- read_envi_header(get_hdr_name(bp_var_sd_path[[parm]]))
        hdr$`band names` <- paste('{',parm,'}',sep = '')
        write_envi_header(hdr, get_hdr_name(bp_var_sd_path[[parm]]))
        bp_var_path[[parm]] <- paste0(bp_var_path[[parm]],'.envi')
        bp_var_sd_path[[parm]] <- paste0(bp_var_sd_path[[parm]],'.envi')
      }
      if (filetype %in% c('GTiff', 'COG')){
        bp_var_path[[parm]] <- paste0(bp_var_path[[parm]],'.tiff')
        bp_var_sd_path[[parm]] <- paste0(bp_var_sd_path[[parm]],'.tiff')
      }
    }
    print('processing completed')
  }
  return(list('mean' = bp_var_path, 'SD' = bp_var_sd_path))
}
