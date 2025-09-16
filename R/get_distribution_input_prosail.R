#' This function generates distribution of biophysical parameters used as input
#' parameters in prosail
#'
#' @param minval list. minimum value for a list of parameters randomly produced
#' @param maxval list. maximum value for a list of parameters randomly produced
#' @param parm_set list. parameters to be set to a given value
#' @param nb_samples numeric. Number of samples to be generated
#' @param type_distrib list. uniform/Gaussian distribution to be applied.
#' @param mean list. mean value for parameters with Gaussian distribution
#' @param sd list. standard deviation for parameters with Gaussian distribution
#'
#' @return input_prosail list. list of nb_samples input parameters for
#' prosail simulation
#' @importFrom stats runif rnorm sd
#' @importFrom truncnorm rtruncnorm
#' @export

get_distribution_input_prosail <- function(minval = NULL, maxval = NULL,
                                           parm_set = NULL, nb_samples = 2000,
                                           type_distrib = NULL, mean = NULL,
                                           sd = NULL){

  if (is.null(type_distrib))
    type_distrib <- data.frame('chl'='Uniform', 'car'='Uniform',
                              'ant' = 'Uniform', 'ewt' = 'Uniform',
                              'lma' = 'Uniform', 'brown'='Uniform',
                              'prot' = 'Uniform', 'cbc' = 'Uniform',
                              'n_struct' = 'Uniform', 'psoil' = 'Uniform',
                              'soil_brightness' = 'Uniform',
                              'lidf_a' = 'Uniform', 'lai' = 'Uniform',
                              'hotspot'='Uniform', 'tto' = 'Uniform',
                              'tts' = 'Uniform', 'psi' = 'Uniform')
  # define all input parameters from PROSAIL
  NAs <- NA*vector(length = nb_samples)
  input_prosail <- data.frame('chl' = NAs, 'car' = NAs, 'ant' = NAs,
                              'brown' = NAs, 'ewt' = NAs, 'lma' = NAs,
                              'prot' = NAs, 'cbc' = NAs, 'n_struct' = NAs,
                              'alpha' = NAs, 'lidf_a' = NAs, 'lidf_b' = NAs,
                              'lai' = NAs, 'hotspot' = NAs, 'tts' = NAs, 'tto' = NAs,
                              'psi' = NAs, 'psoil' = NAs, 'type_lidf' = NAs,
                              'soil_brightness' = NAs)
  input_prosail_names <- names(input_prosail)

  default <- data.frame('chl' = 0, 'car' = 0, 'ant' = 0, 'brown' = 0, 'ewt' = 0,
                        'lma' = 0, 'prot' = 0, 'cbc' = 0, 'n_struct' = 1.5,
                        'alpha' = 40, 'lidf_a' = 0, 'lidf_b' = 0, 'lai' = 2,
                        'hotspot' = 0, 'tts' = 0, 'tto' = 0, 'psi' = 0, 'psoil' = 1,
                        'type_lidf' = 2, 'soil_brightness' = 1)

  # which input parameters should be randomly sampled?
  ParmRand <- input_prosail_names[which(is.element(input_prosail_names,
                                                   names(minval))==TRUE)]
  # which input parameters should be set to a constant value?
  ParmCte <- input_prosail_names[which(is.element(input_prosail_names,
                                                  names(parm_set))==TRUE)]
  # if parameters defined to be sampled randomly and set to constant value
  if (length(intersect(ParmRand,ParmCte))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(intersect(ParmRand,ParmCte))
    message('will be randomly sampled based on minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  all_parms <- c(ParmRand,ParmCte)
  set_to_default <- c()
  if (length(setdiff(input_prosail_names,all_parms))>0){
    # eliminate lma when using PROSPECT-PRO to avoid warnings
    if (!is.na(match('lma', setdiff(input_prosail_names,all_parms)))){
      if (!is.na(match('prot', all_parms)) | !is.na(match('cbc', all_parms))) {
        input_prosail_names <- input_prosail_names[-which(input_prosail_names=='lma')]
        input_prosail$lma <- NULL
      }
    }
    message('WARNING: elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(setdiff(input_prosail_names,all_parms))
    message('These parameters will be set to their default value')
    print(default[setdiff(input_prosail_names,all_parms)])
    set_to_default <- setdiff(input_prosail_names,all_parms)
  }

  # define input_prosail # 1 default value
  if (length(set_to_default)>0){
    for (i in set_to_default)
      input_prosail[[i]] <- default[,i]
  }
  # define input_prosail # 2 Set parameters
  if (length(parm_set)>0){
    for (i in seq_len(length(parm_set))){
      sel <- which(input_prosail_names==names(parm_set)[i])
      input_prosail[[sel]] <- parm_set[,i]
    }
  }
  # define input_prosail # 3 random parameters
  for (i in seq_len(length(minval))){
    sel <- names(minval)[i]
    # if uniform distribution
    if (type_distrib[[sel]] == 'Uniform')
      input_prosail[[sel]] <- runif(nb_samples,
                                    min = minval[1,i],
                                    max = maxval[1,i])
    # if Gaussian distribution
    if (type_distrib[[sel]] == 'Gaussian'){
      set.seed(42)
      input_prosail[[sel]] <- truncnorm::rtruncnorm(n = nb_samples,
                                                    a = minval[[sel]],
                                                    b = maxval[[sel]],
                                                    mean = mean[[sel]],
                                                    sd = sd[[sel]])
    }
  }
  return(input_prosail)
}
