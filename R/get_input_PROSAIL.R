#' generate input_prosail, following
#' - either distribution defined in ATBD
#' - or distribution defined in user
#'
#' @param atbd boolean. TRUE to apply input parameter distribution from ATBD
#' @param GeomAcq list. geom of acquisiton: min and max val for tts, tto & psi
#' @param codistribution_lai boolean. set TYRUE if codistribution with LAI accounted for
#' @param minval list. min val for input parameters
#' @param maxval list. max val for input parameters
#' @param TypeDistrib  list. Type of distribution: 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean & SD for parms sampled with gaussian dist
#' @param parm_set list. list of input parameters set to a specific value
#' @param nb_samples numeric. number of samples in training LUT
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#'
#' @return input_prosail
#' @export

get_input_PROSAIL <- function(atbd = FALSE, GeomAcq = NULL,
                              codistribution_lai = TRUE, minval = NULL,
                              maxval = NULL, TypeDistrib = NULL,
                              GaussianDistrib = NULL, parm_set = NULL,
                              nb_samples = 2000, verbose = FALSE){

  # default parameter values
  defaultVal <- data.frame('CHL' = 40, 'CAR' = 10, 'ANT' = 0, 'EWT' = 0.01,
                           'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5, 'psoil' = 0.5,
                           'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1,
                           'tto' = 0, 'tts' = 30, 'psi' = 80)
  ListParms <- names(defaultVal)

  if (atbd==TRUE){
    #__________________________________________________________________________#
    #                         distribution defined in S2 ATBD                 ##
    #__________________________________________________________________________#
    if (!is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('using PROSAIL parameter distribution defined in S2 ATBD')
        message('http://step.esa.int/docs/extra/ATBD_S2ToolBox_V2.1.pdf')
      }
    }
    input_prosail <- get_atbd_LUT_input(nb_samples = nb_samples,
                                        GeomAcq = GeomAcq,
                                        codistribution_lai = codistribution_lai)
  } else {
    #__________________________________________________________________________#
    #                     user defined range and distribution                 ##
    #__________________________________________________________________________#
    # check consistency between user-defined variables
    # use default distribution if needed
    res <- get_default_LUT_input(TypeDistrib = TypeDistrib,
                                 GaussianDistrib = GaussianDistrib,
                                 minval = minval, maxval = maxval)
    TypeDistrib <- res$TypeDistrib
    GaussianDistrib <- res$GaussianDistrib
    minval <- res$minval
    maxval <- res$maxval

    # fixed parameters
    if (is.na(match('TypeLidf', names(parm_set)))){
      if (is.null(parm_set)){
        parm_set <- data.frame('TypeLidf' = 2)
      } else {
        parm_set <- data.frame(parm_set, 'TypeLidf' = 2)
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
                                                    TypeDistrib = TypeDistrib,
                                                    Mean = GaussianDistrib$Mean,
                                                    Std = GaussianDistrib$Std)
  }
  input_prosail <- data.frame(input_prosail)
  return(input_prosail)
}
