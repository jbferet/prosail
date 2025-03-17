#' generate InputPROSAIL, following
#' - either distribution defined in ATBD
#' - or distribution defined in user
#'
#' @param atbd boolean. TRUE to apply input parameter distribution from ATBD
#' @param GeomAcq list. geom of acquisiton: min and max val for tts, tto & psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#' @param minval list. min val for input parameters
#' @param maxval list. max val for input parameters
#' @param TypeDistrib  list. Type of distribution: 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean & SD for parms sampled with gaussian dist
#' @param ParmSet list. list of input parameters set to a specific value
#' @param nbSamples numeric. number of samples in training LUT
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#'
#' @return InputPROSAIL
#' @export

get_input_PROSAIL <- function(atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                             minval = NULL, maxval = NULL,
                             TypeDistrib = NULL, GaussianDistrib = NULL,
                             ParmSet = NULL, nbSamples = 2000, verbose = FALSE){

  # default parameter values
  ListParms <- c('CHL', 'CAR', 'ANT', 'EWT', 'LMA', 'BROWN', 'N', 'psoil',
                 'LIDFa', 'lai', 'q', 'tto', 'tts', 'psi')
  defaultVal <- data.frame('CHL'=40, 'CAR'=10, 'ANT' = 0, 'EWT' = 0.01,
                           'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5, 'psoil' = 0.5,
                           'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1,
                           'tto' = 0, 'tts' = 30, 'psi' = 80)

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
    InputPROSAIL <- get_atbd_LUT_input(nbSamples = nbSamples,
                                       GeomAcq = GeomAcq,
                                       Codist_LAI = Codist_LAI)
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
    if (is.na(match('TypeLidf', names(ParmSet)))){
      if (is.null(ParmSet)){
        ParmSet <- data.frame('TypeLidf' = 2)
      } else {
        ParmSet <- data.frame(ParmSet, 'TypeLidf' = 2)
      }
    }
    if (is.na(match('alpha', names(ParmSet)))){
      if (is.null(ParmSet)){
        ParmSet <- data.frame('alpha' = 40)
      } else {
        ParmSet <- data.frame(ParmSet, 'alpha' = 40)
      }
    }
    # produce input parameters distribution
    InputPROSAIL <- get_distribution_input_prosail(minval, maxval,
                                                   ParmSet, nbSamples,
                                                   TypeDistrib = TypeDistrib,
                                                   Mean = GaussianDistrib$Mean,
                                                   Std = GaussianDistrib$Std)
  }
  InputPROSAIL <- data.frame(InputPROSAIL)
  return(InputPROSAIL)
}
