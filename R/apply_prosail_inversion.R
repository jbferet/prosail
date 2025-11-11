#' This function applies SVR model on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param mask_path character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param hybrid_model list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param output_dir character. path for directory where results are written
#' @param selected_bands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param band_names character. spectral bands corresponding to the raster
#' @param options
#' - multiplying_factor numeric. to be applied to reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' - tiling boolean. process per block (using raster) or per tile (using terra)
#' - tile_size numeric. size of individual tiles in meters
#' - maxRows numeric. nb of rows to be read from raster
#' - progressBar boolean. if progress bar required
#' - filetype character. driver for raster file
#'
#' @return res character. path for output files corresponding to
#' biophysical properties
#' @importFrom tools file_path_sans_ext
#' @export

apply_prosail_inversion <- function(raster_path, mask_path = NULL, hybrid_model,
                                    output_dir, band_names, selected_bands,
                                    options = NULL){


  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # get raster name
  raster_name <- tools::file_path_sans_ext(basename(raster_path))
  # list of biophysical variables to compute
  bp_vars <- names(hybrid_model)
  options <- set_options_prosail(fun = 'apply_prosail_inversion',
                                 options = options)

  if (options$tiling == FALSE)
    bp_path <- apply_prosail_inversion_per_block(raster_path = raster_path,
                                                 mask_path = mask_path,
                                                 hybrid_model = hybrid_model,
                                                 output_dir = output_dir,
                                                 band_names = band_names,
                                                 selected_bands = selected_bands,
                                                 options = options)

  if (options$tiling == TRUE)
    bp_path <- apply_prosail_inversion_per_tile(raster_path = raster_path,
                                                mask_path = mask_path,
                                                hybrid_model = hybrid_model,
                                                output_dir = output_dir,
                                                band_names = band_names,
                                                selected_bands = selected_bands,
                                                options = options)

  print('processing completed')
  return(bp_path)
}

#' @rdname prosail-deprecated
#' @export
Apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
                                    SelectedBands, bandname, MaskRaster = NULL,
                                    MultiplyingFactor = 10000, maxRows = 100,
                                    bigRaster = FALSE, progressBar = TRUE,
                                    filetype = 'GTiff'){
  .Deprecated("apply_prosail_inversion")
  options <- list('multiplying_factor' = MultiplyingFactor, 'maxRows' = maxRows,
                  'progressBar' = progressBar, 'filetype' = filetype)

  apply_prosail_inversion(raster_path = raster_path,
                          mask_path = MaskRaster,
                          hybrid_model = HybridModel,
                          output_dir = PathOut,
                          band_names = bandname,
                          selected_bands = SelectedBands,
                          options = options)
}
