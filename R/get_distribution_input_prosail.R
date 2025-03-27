#' This function generates distribution of biophysical parameters used as input
#' parameters in PRO4SAIL
#'
#' @param minval list. minimum value for a list of parameters randomly produced
#' @param maxval list. maximum value for a list of parameters randomly produced
#' @param parm_set list. parameters to be set to a given value
#' @param nb_samples numeric. Number of samples to be generated
#' @param TypeDistrib list. uniform/Gaussian distribution to be applied.
#' @param Mean list. mean value for parameters with Gaussian distribution
#' @param Std list. standard deviation for parameters with Gaussian distribution
#'
#' @return input_prosail list. list of nb_samples input parameters for
#' PRO4SAIL simulation
#' @importFrom stats runif rnorm sd
#' @importFrom truncnorm rtruncnorm
#' @export

get_distribution_input_prosail <- function(minval = NULL, maxval = NULL,
                                           parm_set = NULL, nb_samples = 2000,
                                           TypeDistrib = NULL, Mean = NULL,
                                           Std = NULL){

  if (is.null(TypeDistrib))
    TypeDistrib <- data.frame('CHL'='Uniform', 'CAR'='Uniform',
                              'ANT' = 'Uniform', 'EWT' = 'Uniform',
                              'LMA' = 'Uniform', 'BROWN'='Uniform',
                              'PROT' = 'Uniform', 'CBC' = 'Uniform',
                              'N' = 'Uniform', 'psoil' = 'Uniform',
                              'LIDFa' = 'Uniform', 'lai' = 'Uniform',
                              'q'='Uniform', 'tto' = 'Uniform',
                              'tts' = 'Uniform', 'psi' = 'Uniform')
  # define all input parameters from PROSAIL
  fillNA <- NA*vector(length = nb_samples)
  input_prosail <- data.frame('CHL' = fillNA, 'CAR' = fillNA, 'ANT' = fillNA,
                              'BROWN' = fillNA, 'EWT' = fillNA, 'LMA' = fillNA,
                              'PROT' = fillNA, 'CBC' = fillNA, 'N' = fillNA,
                              'alpha' = fillNA, 'LIDFa' = fillNA,
                              'LIDFb' = fillNA, 'lai' = fillNA, 'q' = fillNA,
                              'tts' = fillNA, 'tto' = fillNA, 'psi' = fillNA,
                              'psoil' = fillNA, 'TypeLidf' = fillNA)
  input_prosail_names <- names(input_prosail)

  Default <- data.frame('CHL' = 0, 'CAR' = 0, 'ANT' = 0, 'BROWN' = 0, 'EWT' = 0,
                        'LMA' = 0, 'PROT' = 0, 'CBC' = 0, 'N' = 1.5,
                        'alpha' = 40, 'LIDFa' = 0, 'LIDFb' = 0, 'lai' = 2,
                        'q' = 0, 'tts' = 0, 'tto' = 0, 'psi' = 0, 'psoil' = 1,
                        'TypeLidf' = 2)

  # which input parameters should be randomly sampled?
  ParmRand <- input_prosail_names[which(is.element(input_prosail_names,names(minval))==TRUE)]
  # which input parameters should be set to a constant value?
  ParmCte <- input_prosail_names[which(is.element(input_prosail_names,names(parm_set))==TRUE)]
  # if parameters defined to be sampled randomly and set to constant value
  if (length(intersect(ParmRand,ParmCte))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(intersect(ParmRand,ParmCte))
    message('will be randomly sampled based on the min and max values defined by minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  AllParm <- c(ParmRand,ParmCte)
  Set2Default <- c()
  if (length(setdiff(input_prosail_names,AllParm))>0){
    # eliminate LMA when using PROSPECT-PRO to avoid warnings
    if (!is.na(match('LMA', setdiff(input_prosail_names,AllParm)))){
      if (!is.na(match('PROT', AllParm)) | !is.na(match('CBC', AllParm))) {
        input_prosail_names <- input_prosail_names[-which(input_prosail_names=='LMA')]
        input_prosail$LMA <- NULL
      }
    }
    message('WARNING: elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(setdiff(input_prosail_names,AllParm))
    message('These parameters will be set to their default value')
    print(Default[setdiff(input_prosail_names,AllParm)])
    Set2Default <- setdiff(input_prosail_names,AllParm)
  }

  # define input_prosail # 1 default value
  if (length(Set2Default)>0){
    for (i in Set2Default) input_prosail[[i]] <- Default[,i]
  }
  # define input_prosail # 2 Set parameters
  if (length(parm_set)>0){
    for (i in seq_len(length(parm_set))){
      Sel <- which(input_prosail_names==names(parm_set)[i])
      input_prosail[[Sel]] <- parm_set[,i]
    }
  }
  # define input_prosail # 3 random parameters
  for (i in seq_len(length(minval))){
    Sel <- names(minval)[i]
    # if uniform distribution
    if (TypeDistrib[[Sel]] == 'Uniform')
      input_prosail[[Sel]] <- runif(nb_samples,
                                    min = minval[1,i],
                                    max = maxval[1,i])
    # if Gaussian distribution
    if (TypeDistrib[[Sel]] == 'Gaussian'){
      set.seed(42)
      input_prosail[[Sel]] <- truncnorm::rtruncnorm(n = nb_samples,
                                                    a = minval[[Sel]],
                                                    b = maxval[[Sel]],
                                                    mean = Mean[[Sel]],
                                                    sd = Std[[Sel]])
    }
  }
  return(input_prosail)
}
