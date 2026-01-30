#' computes biophysical properties from a raster subset based on ML regression model
#'
#' @param individual_plot character. path for a raster file
#' @param plot_id character. path for a raster file
#' @param raster_path character. path for a raster file
#' @param bp_vars character. path for a raster file
#' @param hybrid_model list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param output_dir character. path for directory where results are written
#' @param selected_bands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param band_names character. spectral bands corresponding to the raster
#' @param site_name character. spectral bands corresponding to the raster
#' @param mask_path character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param multiplying_factor numeric. multiplying factor used to write
#' reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' @param filetype character. driver for raster file
#' @param p object
#'
#' @return res character. path for output files corresponding to
#' biophysical properties
#' @importFrom matrixStats rowSds
#' @importFrom terra rast crop values writeRaster
#' @importFrom kernlab predict
#' @importFrom stats predict
#' @importFrom tools file_path_sans_ext
#' @export

get_bp_from_tile <- function(individual_plot, plot_id, raster_path, bp_vars,
                             hybrid_model, output_dir, selected_bands, band_names,
                             site_name, mask_path = NULL, multiplying_factor = 10000,
                             filetype = 'GTiff', p = NULL){

  # read image subset
  rast_obj <- terra::rast(raster_path)
  rast_obj <- terra::crop(x = rast_obj, y = individual_plot)
  rast_val <- terra::values(rast_obj)

  # define valid pixels
  sel_pixels <- seq_len(dim(rast_obj)[1]*dim(rast_obj)[2])
  mask_obj <- NULL
  if (!is.null(mask_path)){
    mask_obj <- terra::rast(mask_path)
    mask_obj <- terra::crop(x = mask_obj, y = individual_plot)
    sel_pixels <- which(terra::values(mask_obj)==1)
  }
  rast_val <- rast_val[sel_pixels, ]

  # define path for maps produced for each biophysical variable
  bp_var_path <- bp_var_sd_path <- list()
  # for each parameter
  for (parm in bp_vars){
    # define file name and check if exists
    output_dir_mean <- file.path(output_dir, 'mean')
    dir.create(path = output_dir_mean, showWarnings = FALSE, recursive = TRUE)
    bp_var_path[[parm]] <- file.path(output_dir_mean,
                                     paste(site_name, plot_id, parm, sep = '_'))

    output_dir_sd <- file.path(output_dir, 'sd')
    dir.create(path = output_dir_sd, showWarnings = FALSE, recursive = TRUE)
    bp_var_sd_path[[parm]] <- file.path(output_dir_sd, paste(site_name, plot_id,
                                                             parm, 'sd', sep = '_'))
    if (filetype %in% c('GTiff', 'COG')){
      bp_var_path[[parm]] <- paste0(bp_var_path[[parm]],'.tiff')
      bp_var_sd_path[[parm]] <- paste0(bp_var_sd_path[[parm]],'.tiff')
    }

    # if files do not exist
    if (! file.exists(bp_var_path[[parm]]) & ! file.exists(bp_var_sd_path[[parm]])){
      # define spectral bands
      sel_bands <- match(selected_bands[[parm]], band_names)
      block_val <- rast_val[, sel_bands]
      if (!inherits(hybrid_model[[parm]][[1]], what = 'liquidSVM') &
          !inherits(hybrid_model[[parm]][[1]], what = 'ksvm'))
        colnames(block_val) <- colnames(hybrid_model[[parm]][[1]]$trainingData)[-1]

      rast_mean <- rast_sd <- NA*rast_obj[[1]]
      names(rast_mean) <- paste('mean', parm)
      names(rast_sd) <- paste('sd', parm)
      if (length(sel_pixels)>0){
        block_val <- block_val/multiplying_factor
        estimates <- list()
        for (modind in seq_len(length(hybrid_model[[parm]]))){
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
        rast_mean[sel_pixels] <- mean_estimate
        rast_sd[sel_pixels] <- sd_estimate
      }

      terra::writeRaster(x = rast_mean, filename = bp_var_path[[parm]],
                         filetype = filetype, overwrite = TRUE,
                         gdal=c("COMPRESS=DEFLATE", "TFW=YES"))
      terra::writeRaster(x = rast_sd, filename = bp_var_sd_path[[parm]],
                         filetype = filetype, overwrite = TRUE,
                         gdal=c("COMPRESS=DEFLATE", "TFW=YES"))

      # write biophysical variable name in headers if ENVI type
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
    }
  }
  if (!is.null(p))
    p()
  return(list('mean' = bp_var_path, 'sd' = bp_var_sd_path))
}
