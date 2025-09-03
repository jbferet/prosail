#' generate input_prosail, following
#' - either distribution defined in ATBD
#' - or distribution defined in user
#'
#' @param atbd boolean. TRUE to apply input parameter distribution from ATBD
#' @param geom_acq list. geom of acquisiton: min and max val for tts, tto & psi
#' @param codistribution_lai boolean. TRUE accounts for codistribution with LAI
#' @param minval list. min val for input parameters
#' @param maxval list. max val for input parameters
#' @param type_distrib  list. Type of distribution: 'Uniform' or 'Gaussian'
#' @param gaussian_distrib  list. Mean & SD for parms sampled with gaussian dist
#' @param parm_set list. list of input parameters set to a specific value
#' @param nb_samples numeric. number of samples in training LUT
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#'
#' @return input_prosail
#' @export

get_input_prosail <- function(atbd = FALSE, geom_acq = NULL,
                              codistribution_lai = TRUE, minval = NULL,
                              maxval = NULL, type_distrib = NULL,
                              gaussian_distrib = NULL, parm_set = NULL,
                              nb_samples = 2000, verbose = FALSE){

  if (atbd==TRUE | tolower(x = atbd)=='v2'){
    #__________________________________________________________________________#
    #                         distribution defined in S2 ATBD                 ##
    #__________________________________________________________________________#
    if (!is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('using PROSAIL parameter distribution defined in S2 ATBD v2')
        message('http://step.esa.int/docs/extra/ATBD_S2ToolBox_V2.1.pdf')
      }
    }
    input_prosail <- get_atbd_lut_input(nb_samples = nb_samples,
                                        geom_acq = geom_acq,
                                        codistribution_lai = codistribution_lai)
  } else if (tolower(x = atbd)=='v3'){
    #__________________________________________________________________________#
    #                         distribution defined in S2 ATBD                 ##
    #__________________________________________________________________________#
    input_prosail <- get_atbd_v3_lut_input(nb_samples = nb_samples,
                                           geom_acq = geom_acq,
                                           codistribution_lai = codistribution_lai)
  } else {
    #__________________________________________________________________________#
    #                     user defined range and distribution                 ##
    #__________________________________________________________________________#
    # check consistency between user-defined variables
    # use default distribution if needed
    res <- get_default_lut_input(type_distrib = type_distrib,
                                 gaussian_distrib = gaussian_distrib,
                                 minval = minval, maxval = maxval)
    type_distrib <- res$type_distrib
    gaussian_distrib <- res$gaussian_distrib
    minval <- res$minval
    maxval <- res$maxval

    # fixed parameters
    if (is.na(match('type_lidf', names(parm_set)))){
      if (is.null(parm_set)){
        parm_set <- data.frame('type_lidf' = 2)
      } else {
        parm_set <- data.frame(parm_set, 'type_lidf' = 2)
      }
    }
    if (is.na(match('alpha', names(parm_set)))){
      if (is.null(parm_set)){
        parm_set <- data.frame('alpha' = 40)
      } else {
        parm_set <- data.frame(parm_set, 'alpha' = 40)
      }
    }
    # produce input parameters distribution
    input_prosail <- get_distribution_input_prosail(minval, maxval,
                                                    parm_set, nb_samples,
                                                    type_distrib = type_distrib,
                                                    mean = gaussian_distrib$mean,
                                                    sd = gaussian_distrib$sd)
  }
  input_prosail <- data.frame(input_prosail)
  return(input_prosail)
}


#' @rdname prosail-deprecated
#' @export
get_InputPROSAIL <- function(atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                             minval = NULL, maxval = NULL,
                             TypeDistrib = NULL, GaussianDistrib = NULL,
                             ParmSet = NULL, nbSamples = 2000, verbose = FALSE){
  .Deprecated("get_input_prosail")
  get_input_prosail(atbd, GeomAcq, Codist_LAI, minval, maxval, TypeDistrib,
                    GaussianDistrib, ParmSet, nbSamples, verbose)
}
