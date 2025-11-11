#' set options
#'
#' @param fun character. name of the function which has optional parameters
#' @param options list. including
#' - cloudcover
#' - path_S2tilinggrid
#' - overwrite
#' - geomAcq
#' - additional_process
#' - bands2correct
#' - fraction_vegetation
#' - RadiometricFilter
#'
#' @return options with default values when missing
#' @export

set_options_prosail <- function(fun, options = NULL){

  if (fun == 'train_prosail_inversion'){
    if (is.null(options$codistribution_lai))
      options$codistribution_lai <- TRUE
    if (is.null(options$type_distrib))
      options$type_distrib <- NULL
    if (is.null(options$gaussian_distrib))
      options$gaussian_distrib <- NULL
    if (is.null(options$parm_set))
      options$parm_set <- NULL
    if (is.null(options$SAILversion))
      options$SAILversion <- '4SAIL'
    if (is.null(options$brown_lop))
      options$brown_lop <- NULL
    if (is.null(options$nb_samples))
      options$nb_samples <- 2000
    if (is.null(options$nb_models))
      options$nb_models <- 20
    if (is.null(options$replacement))
      options$replacement <- TRUE
    if (is.null(options$noise_level))
      options$noise_level <- NULL
    if (is.null(options$spec_prospect))
      options$spec_prospect <- NULL
    if (is.null(options$spec_soil))
      options$spec_soil <- NULL
    if (is.null(options$spec_atm))
      options$spec_atm <- NULL
    if (is.null(options$method))
      options$method <- 'liquidSVM'
    if (is.null(options$verbose))
      options$verbose <- FALSE
  }
  if (fun == 'apply_prosail_inversion'){
    if (is.null(options$multiplying_factor))
      options$multiplying_factor <- 10000
    if (is.null(options$maxRows))
      options$maxRows <- 100
    if (is.null(options$progressBar))
      options$progressBar <- TRUE
    if (is.null(options$filetype))
      options$filetype <- 'GTiff'
    if (is.null(options$tiling))
      options$tiling <- FALSE
    if (is.null(options$tile_size))
      options$tile_size <- 5000
  }
  if (fun == 'apply_prosail_inversion_per_tile'){
    if (is.null(options$multiplying_factor))
      options$multiplying_factor <- 10000
    if (is.null(options$filetype))
      options$filetype <- 'GTiff'
    if (is.null(options$tiling))
      options$tiling <- FALSE
    if (is.null(options$tile_size))
      options$tile_size <- 5000
  }
  if (fun == 'apply_prosail_inversion_per_block'){
    if (is.null(options$multiplying_factor))
      options$multiplying_factor <- 10000
    if (is.null(options$maxRows))
      options$maxRows <- 100
    if (is.null(options$progressBar))
      options$progressBar <- TRUE
    if (is.null(options$filetype))
      options$filetype <- 'GTiff'

  }
  return(options)
}

