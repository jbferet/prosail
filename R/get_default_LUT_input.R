#' Sets default values for PROSAIL LUT simulation when not defined by user
#'
#' @param TypeDistrib list. uniform/Gaussian distribution to apply.
#' @param GaussianDistrib  list. Mean & STD for parms with gaussian distribution
#' @param minval list. minimum value for a list of parameters randomly produced
#' @param maxval list. maximum value for a list of parameters randomly produced
#'
#' @return res list. default values corresponding to NULL input parameters
#' @export

get_default_LUT_input <- function(TypeDistrib = NULL,
                                  GaussianDistrib = NULL,
                                  minval = NULL,
                                  maxval = NULL){

  # define mean / sd for gaussian
  if (is.null(GaussianDistrib)) GaussianDistrib <- list('Mean'=NULL,'Std'=NULL)
  # check consistency between TypeDistrib and GaussianDistrib
  if (!is.null(TypeDistrib)){
    namesDist <- names(TypeDistrib)
    whichGauss <- which(TypeDistrib=='Gaussian')
    if (length(whichGauss)>0){
      namesGauss <- names(GaussianDistrib$Mean)
      selGauss <- namesDist[whichGauss]
      matchGauss <- match(selGauss, namesGauss)
      WhichMissed <- selGauss[which(is.na(matchGauss))]
      if (length(WhichMissed)>0){
        message(paste('missing Mean and Std for GaussianDistrib of parameter',
                      WhichMissed))
        stop(message('abort process'))
      }
    }
  }

  namesParms <- names(TypeDistrib)
  namesMin <- names(minval)
  namesMax <- names(maxval)
  MatchingVars <- c(match(namesParms, namesMin), match(namesMin, namesParms),
                    match(namesParms, namesMax), match(namesMax, namesParms),
                    match(namesMin, namesMax), match(namesMax, namesMin))
  if (length(which(is.na(MatchingVars)))>0){
    message('Make sure TypeDistrib, minval & maxval share the same parameters')
    stop(message('abort process'))
  }
  # define uniform / gaussian distribution
  if (is.null(TypeDistrib))
    TypeDistrib <- data.frame('CHL'='Uniform', 'CAR'='Uniform', 'ANT'='Uniform',
                              'BROWN'='Uniform', 'EWT' = 'Uniform',
                              'LMA' = 'Uniform', 'N' = 'Uniform',
                              'psoil' = 'Uniform', 'LIDFa' = 'Uniform',
                              'lai' = 'Uniform', 'q'='Uniform',
                              'tto' = 'Uniform','tts' = 'Uniform',
                              'psi' = 'Uniform')
  # define min and max values
  if (is.null(minval))
    minval <- data.frame('CHL' = 10, 'CAR' = 0, 'EWT' = 0.01, 'ANT' = 0,
                         'LMA' = 0.005, 'N' = 1.0, 'psoil' = 0.0, 'BROWN'=0.0,
                         'LIDFa' = 20, 'lai' = 0.5, 'q'=0.1, 'tto' = 0,
                         'tts' = 20, 'psi' = 80)
  if (is.null(maxval))
    maxval <- data.frame('CHL' = 75, 'CAR' = 15, 'EWT' = 0.03, 'ANT' = 2,
                         'LMA' = 0.03, 'N' = 2.0, 'psoil' = 1.0, 'BROWN'=0.5,
                         'LIDFa' = 70, 'lai' = 7, 'q'=0.2, 'tto' = 5,
                         'tts' = 30, 'psi' = 110)
  res <- list('TypeDistrib' = TypeDistrib,
              'GaussianDistrib' = GaussianDistrib,
              'minval' = minval, 'maxval' = maxval)
  return(res)
}
