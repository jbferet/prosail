#' This function applies ML regression on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param mask_path character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param hybrid_model list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param output_dir character. path for directory where results are written
#' @param band_names character. spectral bands corresponding to the raster
#' @param selected_bands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param options
#' - multiplying_factor numeric. to be applied to reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' - maxRows numeric. nb of rows to be read from raster
#' - progressBar boolean. if progress bar required
#' - filetype character. driver for raster file
#'
#' @return res character. path for output files corresponding to
#' biophysical properties
#' @importFrom progressr handlers progressor with_progress
#' @importFrom raster raster brick blockSize readStart readStop getValues
#' writeStart writeStop writeValues
#' @importFrom matrixStats rowSds
#' @importFrom tools file_path_sans_ext
#' @importFrom utils installed.packages
#' @export
#'
apply_prosail_inversion_per_block <- function(raster_path, mask_path = NULL,
                                              hybrid_model, output_dir,
                                              band_names, selected_bands,
                                              options = NULL){
  options <- set_options_prosail(fun = 'apply_prosail_inversion_per_block',
                                 options = options)
  multiplying_factor <- options$multiplying_factor
  maxRows<- options$maxRows
  progressBar <- options$progressBar
  filetype <- options$filetype

  bp_vars <- names(hybrid_model)
  raster_name <- tools::file_path_sans_ext(basename(raster_path))
  # define path for maps produced for each biophysical variable
  bp_var_path <- bp_var_sd_path <- list()
  # check if bigRaster supported
  for (parm in bp_vars){
    print(paste('Computing',parm,sep = ' '))
    # which bands to select
    sel_bands <- match(selected_bands[[parm]], band_names)
    # read by chunk to avoid memory problem
    blk <- raster::blockSize(raster::brick(raster_path))
    if (blk$nrows[1]>maxRows){
      nrows_indiv <- maxRows
      blkud <- list()
      blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
      if (sum(blkud$nrows)<sum(blk$nrows))
        blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
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
    if (!is.null(mask_path)){
      if (file.exists(mask_path)){
        r_inmask <- raster::readStart(raster(mask_path))
      } else if (!file.exists(mask_path)){
        message('WARNING: Mask file does not exist:')
        print(mask_path)
      }
    }
    if (progressBar == TRUE){
      # initiate progress bar
      pgbarlength <- length(hybrid_model[[parm]])*blk$n
      pb <- progress_bar$new(
        format = "Mapping from raster [:bar] :percent in :elapsedfull , estimated time remaining :eta",
        total = pgbarlength, clear = FALSE, width= 100)
    }
    # output files
    formatfile <- filetype
    if (filetype == 'EHdr')
      formatfile <- 'ENVI'
    bp_var_path[[parm]] <- file.path(output_dir, paste(raster_name, parm,
                                                       sep = '_'))
    bp_var_sd_path[[parm]] <- file.path(output_dir,paste(raster_name, parm,
                                                         'STD', sep = '_'))

    if (filetype %in% c('GTiff', 'COG')){
      bp_var_path[[parm]] <- paste0(bp_var_path[[parm]],'.tiff')
      bp_var_sd_path[[parm]] <- paste0(bp_var_sd_path[[parm]],'.tiff')
    }

    r_out_mean <- writeStart(raster(raster_path),
                             filename = bp_var_path[[parm]],
                             format = formatfile, overwrite = TRUE)
    r_out_sd <- writeStart(raster(raster_path),
                           filename = bp_var_sd_path[[parm]],
                           format = formatfile, overwrite = TRUE)

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
        MaskVal <- raster::getValues(r_inmask, row = blk$row[i],
                                     nrows = blk$nrows[i])
        select_pixels <- which(MaskVal ==1)
        block_val <- block_val[select_pixels,sel_bands]
      }
      # add name for variables

      if (!inherits(hybrid_model[[parm]][[1]], what = 'liquidSVM') &
          !inherits(hybrid_model[[parm]][[1]], what = 'ksvm'))
        colnames(block_val) <- colnames(hybrid_model[[parm]][[1]]$trainingData)[-1]

      mean_estimate_full <- NA*vector(length = full_length)
      sd_estimate_full <- NA*vector(length = full_length)
      if (length(select_pixels)>0){
        block_val <- block_val/multiplying_factor
        estimates <- list()
        for (modind in seq_len(length(hybrid_model[[parm]]))){
          if (progressBar == TRUE)
            pb$tick()
          if (!inherits(hybrid_model[[parm]][[1]], what = 'ksvm')){
            estimates[[modind]] <- stats::predict(hybrid_model[[parm]][[modind]],
                                                  block_val)
          } else if (inherits(hybrid_model[[parm]][[1]], what = 'ksvm')){
            estimates[[modind]] <- kernlab::predict(object = hybrid_model[[parm]][[modind]],
                                                    block_val)
          }
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
    print('processing completed')
  }
  return(list('mean' = bp_var_path, 'SD' = bp_var_sd_path))
}
