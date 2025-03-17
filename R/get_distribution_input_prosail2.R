#' This function generates distribution of biophysical parameters used as input
#' parameters in PRO4SAIL
#'
#' @param minval list. min val set for list of parameters randomly produced
#' @param maxval list. max val set for list of parameters randomly produced
#' @param ParmSet list.Defines the parameters to be set to a given value
#' @param nbSamples numeric. Number of samples to be generated
#' @param TypeDistrib list. uniform / Gaussian distribution to be applied
#' @param Mean list. mean value for parameters with Gaussian distribution
#' @param Std list. standard deviation for parameters with Gaussian distribution
#'
#' @return InputPROSAIL list. nbSamples input parameters for PRO4SAIL simulation
#' @importFrom stats runif rnorm sd
#' @importFrom truncnorm rtruncnorm
#' @export
get_distribution_input_prosail2 <- function(minval,maxval,ParmSet,nbSamples,
                                            TypeDistrib = NULL, Mean = NULL,
                                            Std = NULL){

  if (is.null(TypeDistrib)){
    TypeDistrib <- data.frame('CHL'='Uniform', 'CAR'='Uniform',
                              'ANT' = 'Uniform', 'EWT' = 'Uniform',
                              'LMA' = 'Uniform', 'BROWN'='Uniform',
                              'PROT' = 'Uniform', 'CBC' = 'Uniform',
                              'N' = 'Uniform', 'psoil' = 'Uniform',
                              'LIDFa' = 'Uniform', 'lai' = 'Uniform',
                              'q'='Uniform', 'tto' = 'Uniform',
                              'tts' = 'Uniform', 'psi' = 'Uniform',
                              'fraction_brown' = 'Uniform', 'diss' = 'Uniform',
                              'Cv' = 'Uniform', 'Zeta' = 'Uniform')
  }
  # define all input parameters from PROSAIL
  InVar <- c('CHL','CAR','ANT','BROWN','EWT','LMA',
             'PROT','CBC','N','alpha','LIDFa','LIDFb',
             'lai','q','tts','tto','psi','psoil','TypeLidf',
             'fraction_brown','diss','Cv','Zeta')

  InputPROSAIL <- list('CHL' = c(), 'CAR' = c(), 'ANT' = c(), 'BROWN' = c(),
                       'EWT' = c(), 'LMA' = c(), 'PROT' = c(), 'CBC' = c(),
                       'N' = c(), 'alpha' = c(), 'LIDFa' = c(), 'LIDFb' = c(),
                       'lai' = c(), 'q' = c(), 'tts' = c(), 'tto' = c(),
                       'psi' = c(), 'psoil' = c(), 'TypeLidf' = c(),
                       'fraction_brown' = c(), 'diss' = c(), 'Cv' = c(),
                       'Zeta' = c())
  Default <- data.frame('CHL' = 0, 'CAR' = 0, 'ANT'=0, 'BROWN' = 0, 'EWT' = 0,
                        'LMA' = 0, 'PROT' = 0, 'CBC' = 0, 'N' = 1.5,
                        'alpha' = 40, 'LIDFa' = 0, 'LIDFb' = 0, 'lai' = 2,
                        'q' = 0, 'tts' = 0, 'tto' = 0,'psi'=0, 'psoil' = 1,
                        'TypeLidf' = 2, 'fraction_brown' = 0.5, 'diss' = 0.5,
                        'Cv' = 1, 'Zeta' = 1)

  # which input parameters should be randomly sampled?
  ParmRand <- which(is.element(InVar,names(minval))==TRUE)
  # which input parameters should be set to a constant value?
  ParmCte <- which(is.element(InVar,names(ParmSet))==TRUE)
  # if parameters are defined to be sampled randomly and set to constant value
  if (length(intersect(ParmRand,ParmCte))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(InVar[intersect(ParmRand,ParmCte)])
    message('will be randomly sampled based on the min and max values defined by minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  AllParm <- c(ParmRand,ParmCte)
  FullList <- seq_len(length(InVar))
  Set2Default <- c()
  if (length(setdiff(FullList,AllParm))>0){
    message('WARNING: some elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(InVar[setdiff(FullList,AllParm)])
    message('These parameters will be set to their default value')
    print(Default[setdiff(FullList,AllParm)])
    Set2Default <- setdiff(FullList,AllParm)
  }

  # define InputPROSAIL # 1 default value
  if (length(Set2Default)>0){
    for (i in Set2Default) InputPROSAIL[[i]] <- Default[,i]
  }

  # define InputPROSAIL # 2 Set parameters
  if (length(ParmSet)>0){
    for (i in seq_len(length(ParmSet))){
      Sel <- which(InVar==names(ParmSet)[i])
      InputPROSAIL[[Sel]] <- ParmSet[,i]
    }
  }

  # define InputPROSAIL # 3 random parameters
  for (i in seq_len(length(minval))){
    Sel <- names(minval)[i]
    # if uniform distribution
    if (TypeDistrib[[Sel]] == 'Uniform'){
      InputPROSAIL[[Sel]] <- array(runif(nbSamples,
                                         min = minval[1,i],
                                         max = maxval[1,i]),
                                   dim = c(nbSamples,1))
      # if Gaussian distribution
    } else if (TypeDistrib[[Sel]] == 'Gaussian'){
      set.seed(42)
      InputPROSAIL[[Sel]] <- truncnorm::rtruncnorm(n = nbSamples,
                                                   a = minval[[Sel]],
                                                   b = maxval[[Sel]],
                                                   mean = Mean[[Sel]],
                                                   sd = Std[[Sel]])
      InputPROSAIL[[Sel]] <- array(InputPROSAIL[[Sel]],dim = c(nbSamples,1))
    }
  }
  return(InputPROSAIL)
}
