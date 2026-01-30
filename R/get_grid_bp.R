#' return vector grid corresponding to an area of interest
#'
#' @param raster_path character. path for vector file defining aoi
#' @param cellsize numeric. square cell size in meters
#'
#' @return dsn_grid path for grid of aoi
#' @importFrom terra rast buffer as.polygons
#' @importFrom sf st_as_sf st_multipolygon st_make_grid st_intersects st_buffer sf_use_s2
#' @export

get_grid_bp <- function(raster_path, cellsize = 10000){

  rast_obj <- terra::rast(raster_path)[[1]]
  mask_buffer <- terra::buffer(x = rast_obj , 20)
  mask_polygon <- terra::as.polygons(mask_buffer)
  aoi <- sf::st_as_sf(mask_polygon, quiet = T)
  crs_init <- sf::st_crs(aoi)
  # aoi <- sf::st_multipolygon(lapply(aoi$geometry[2], function(x) x[1]))
  aoi <- sf::st_multipolygon(lapply(aoi$geometry, function(x) x[1]))
  aoi <- sf::st_sf(sf::st_sfc(aoi))
  sf::st_crs(aoi) <- crs_init
  aoi_grid <- aoi |> sf::st_make_grid(cellsize = cellsize, square = TRUE)

  # select grid cells intersecting with initial aoi
  suppressMessages(sf::sf_use_s2(FALSE))
  intersect <- as.data.frame(sf::st_intersects(x = aoi_grid, aoi))
  suppressMessages(sf::sf_use_s2(TRUE))
  aoi_grid <- aoi_grid[intersect$row.id]
  # add buffer to avoid border effects
  aoi_grid <- sf::st_buffer(x = aoi_grid,
                            dist = 10*terra::res(x = rast_obj))
  sf::st_crs(aoi_grid) <- crs_init

  # save grid
  filename <- file.path(dirname(raster_path),
                        paste0('aoi_tiles_', cellsize, 'm.gpkg'))
  sf::st_write(obj = aoi_grid, dsn = filename,
               overwrite = T, append = F,
               driver = 'GPKG', quiet = T)
  nb_tiles <- length(aoi_grid)
  plots <- get_plot_list(dsn = filename, nbdigits = nchar(nb_tiles))
  return(plots)
}
