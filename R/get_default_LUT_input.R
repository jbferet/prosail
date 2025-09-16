#' Sets default values for PROSAIL LUT simulation when not defined by user
#'
#' @param type_distrib list. uniform/Gaussian distribution to apply.
#' @param gaussian_distrib  list. mean & sd for parms with gaussian dist
#' @param minval list. minimum value for a list of parameters randomly produced
#' @param maxval list. maximum value for a list of parameters randomly produced
#'
#' @return res list. default values corresponding to NULL input parameters
#' @export

get_default_lut_input <- function(type_distrib = NULL,
                                  gaussian_distrib = NULL,
                                  minval = NULL,
                                  maxval = NULL){

  # define mean / sd for gaussian
  if (is.null(gaussian_distrib)) gaussian_distrib <- list('mean'=NULL,'sd'=NULL)
  # check consistency between type_distrib and gaussian_distrib
  if (!is.null(type_distrib)){
    names_dist <- names(type_distrib)
    which_gauss <- which(type_distrib=='Gaussian')
    if (length(which_gauss)>0){
      names_gauss <- names(gaussian_distrib$Mean)
      sel_gauss <- names_dist[which_gauss]
      match_gauss <- match(sel_gauss, names_gauss)
      which_missed <- sel_gauss[which(is.na(match_gauss))]
      if (length(which_missed)>0){
        message(paste('missing mean and sd for gaussian_distrib of parameter',
                      which_missed))
        stop(message('abort process'))
      }
    }
  }

  names_parms <- names(type_distrib)
  names_min <- names(minval)
  names_max <- names(maxval)
  matching_vars <- c(match(names_parms, names_min),
                     match(names_min, names_parms),
                     match(names_parms, names_max),
                     match(names_max, names_parms),
                     match(names_min, names_max),
                     match(names_max, names_min))
  if (length(which(is.na(matching_vars)))>0){
    message('Make sure type_distrib, minval & maxval share the same parameters')
    stop(message('abort process'))
  }
  # define uniform / gaussian distribution
  if (is.null(type_distrib))
    type_distrib <- data.frame('chl'='Uniform', 'car'='Uniform', 'ant'='Uniform',
                              'brown'='Uniform', 'ewt' = 'Uniform',
                              'lma' = 'Uniform', 'n_struct' = 'Uniform',
                              'psoil' = 'Uniform', 'lidf_a' = 'Uniform',
                              'lai' = 'Uniform', 'hotspot'='Uniform',
                              'tto' = 'Uniform','tts' = 'Uniform',
                              'psi' = 'Uniform', 'soil_brightness' = 'Uniform')
  # define min and max values
  if (is.null(minval))
    minval <- data.frame('chl' = 10, 'car' = 0, 'ewt' = 0.01, 'ant' = 0,
                         'lma' = 0.005, 'brown'=0.0, 'n_struct' = 1.0,
                         'psoil' = 0.0, 'soil_brightness' = 0.5,
                         'lidf_a' = 20, 'lai' = 0.5, 'hotspot'=0.1,
                         'tto' = 0, 'tts' = 20, 'psi' = 80)
  if (is.null(maxval))
    maxval <- data.frame('chl' = 75, 'car' = 15, 'ewt' = 0.03, 'ant' = 2,
                         'lma' = 0.03, 'brown'=0.5, 'n_struct' = 2.0,
                         'psoil' = 1.0, 'soil_brightness' = 3.5,
                         'lidf_a' = 70, 'lai' = 7, 'hotspot'=0.2,
                         'tto' = 5, 'tts' = 30, 'psi' = 110)
  res <- list('type_distrib' = type_distrib,
              'gaussian_distrib' = gaussian_distrib,
              'minval' = minval, 'maxval' = maxval)
  return(res)
}




#' @rdname prosail-deprecated
#' @export
get_default_LUT_input <- function(TypeDistrib = NULL,
                                  GaussianDistrib = NULL,
                                  minval = NULL,
                                  maxval = NULL){
  .Deprecated("get_default_lut_input")
  get_default_lut_input(TypeDistrib, GaussianDistrib, minval, maxval)
}
