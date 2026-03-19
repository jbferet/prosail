#' This function applies SVR model on raster data in order to estimate
#' vegetation biophysical properties
#'
#' @param raster_path character. path for a raster file
#' @param hybrid_model list. hybrid models produced from train_prosail_inversion
#' each element of the list corresponds to a set of hybrid models for a
#' vegetation parameter
#' @param output_dir character. path for directory where results are written
#' @param selected_bands list. list of spectral bands to be selected from raster
#' (identified by name of vegetation parameter)
#' @param band_names character. spectral bands corresponding to the raster
#' @param mask_path character. path for binary mask defining ON (1) and OFF (0)
#' pixels in the raster
#' @param multiplying_factor numeric. multiplying factor used to write
#' reflectance in the raster
#' --> PROSAIL simulates reflectance between 0 and 1, and raster data expected
#' in the same range
#' @param tile_size numeric. size of individual tiles in meters
#' @param filetype character. driver for raster file
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

apply_prosail_inversion_tiles <- function(raster_path, hybrid_model, output_dir,
                                          selected_bands, band_names, mask_path = NULL,
                                          multiplying_factor = 10000, tile_size = 10000,
                                          filetype = 'GTiff'){

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # get raster name
  raster_name <- tools::file_path_sans_ext(basename(raster_path))
  # list of biophysical variables to compute
  bp_vars <- names(hybrid_model)
  # define grid
  plots <- get_grid_bp(raster_path, cellsize = tile_size)
  site_name <- tools::file_path_sans_ext(basename(raster_name))

  # for each plot
  # nbCPU <- 0
  # if (nbCPU == 0){
  #   cl <- parallel::makeCluster(nbCPU)
  #   plan("cluster", workers = cl)
  #   bp_var_path <- future.apply::future_mapply(FUN = get_bp_from_tile,
  #                                              individual_plot = plots,
  #                                              plot_id = as.list(names(plots)),
  #                                              MoreArgs = list(raster_path = raster_path,
  #                                                              mask_path = mask_path,
  #                                                              bp_vars = bp_vars,
  #                                                              hybrid_model = hybrid_model,
  #                                                              output_dir = output_dir,
  #                                                              selected_bands = selected_bands,
  #                                                              band_names   = band_names,
  #                                                              site_name = site_name,
  #                                                              multiplying_factor = multiplying_factor,
  #                                                              filetype = filetype),
  #                                              future.seed = TRUE)
  #   parallel::stopCluster(cl)
  #   plan(sequential)
  # } else {
    handlers("cli")
    suppressWarnings(with_progress({
      p <- progressr::progressor(steps = length(plots),
                                 message = 'compute BP from tiles')
      bp_var_path <- mapply(FUN = get_bp_from_tile,
                            individual_plot = plots,
                            plot_id = as.list(names(plots)),
                            MoreArgs = list(raster_path = raster_path,
                                            mask_path = mask_path,
                                            bp_vars = bp_vars,
                                            hybrid_model = hybrid_model,
                                            output_dir = output_dir,
                                            selected_bands = selected_bands,
                                            band_names = band_names,
                                            site_name = site_name,
                                            multiplying_factor = multiplying_factor,
                                            filetype = filetype,
                                            p = p),
                            SIMPLIFY = F)}))
  # }
  return(bp_var_path)
}
