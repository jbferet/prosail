#' This function generates distribution of biophysical parameters used as input
#' parameters in prosail
#'
#' @param minval list. min val set for list of parameters randomly produced
#' @param maxval list. max val set for list of parameters randomly produced
#' @param parm_set list.Defines the parameters to be set to a given value
#' @param nb_samples numeric. Number of samples to be generated
#' @param type_distrib list. uniform / Gaussian distribution to be applied
#' @param mean list. mean value for parameters with Gaussian distribution
#' @param sd list. standard deviation for parameters with Gaussian distribution
#'
#' @return input_prosail list. input parameters for prosail simulation
#' @importFrom stats runif rnorm sd
#' @importFrom truncnorm rtruncnorm
#' @export
get_distribution_input_prosail2 <- function(minval, maxval, parm_set,
                                            nb_samples, type_distrib = NULL,
                                            mean = NULL, sd = NULL){

  if (is.null(type_distrib)){
    type_distrib <- data.frame('chl'='Uniform', 'car'='Uniform',
                              'ant' = 'Uniform', 'ewt' = 'Uniform',
                              'lma' = 'Uniform', 'brown'='Uniform',
                              'prot' = 'Uniform', 'cbc' = 'Uniform',
                              'n_struct' = 'Uniform', 'psoil' = 'Uniform',
                              'lidf_a' = 'Uniform', 'lai' = 'Uniform',
                              'q'='Uniform', 'tto' = 'Uniform',
                              'tts' = 'Uniform', 'psi' = 'Uniform',
                              'fraction_brown' = 'Uniform', 'diss' = 'Uniform',
                              'cv' = 'Uniform', 'zeta' = 'Uniform')
  }
  # define all input parameters from PROSAIL
  input_prosail <- list('chl' = c(), 'car' = c(), 'ant' = c(), 'brown' = c(),
                        'ewt' = c(), 'lma' = c(), 'prot' = c(), 'cbc' = c(),
                        'n_struct' = c(), 'alpha' = c(), 'lidf_a' = c(),
                        'lidf_b' = c(), 'lai' = c(), 'q' = c(), 'tts' = c(),
                        'tto' = c(), 'psi' = c(), 'psoil' = c(),
                        'type_lidf' = c(), 'fraction_brown' = c(), 'diss' = c(),
                        'cv' = c(), 'zeta' = c())
  input_prosail_names <- names(input_prosail)

  default <- data.frame('chl' = 0, 'car' = 0, 'ant'=0, 'brown' = 0, 'ewt' = 0,
                        'lma' = 0, 'prot' = 0, 'cbc' = 0, 'n_struct' = 1.5,
                        'alpha' = 40, 'lidf_a' = 0, 'lidf_b' = 0, 'lai' = 2,
                        'q' = 0, 'tts' = 0, 'tto' = 0,'psi'=0, 'psoil' = 1,
                        'type_lidf' = 2, 'fraction_brown' = 0.5, 'diss' = 0.5,
                        'cv' = 1, 'zeta' = 1)

  # which input parameters should be randomly sampled?
  parm_random <- which(is.element(input_prosail_names,names(minval))==TRUE)
  # which input parameters should be set to a constant value?
  parm_constante <- which(is.element(input_prosail_names,names(parm_set))==TRUE)
  # if parameters are defined to be sampled randomly and set to constant value
  if (length(intersect(parm_random,parm_constante))>0){
    message('WARNING: some elements are defined as parameters to sample')
    message('and parameters to set: ')
    print(input_prosail_names[intersect(parm_random,parm_constante)])
    message('will be randomly sampled based on minval and maxval')
  }
  # if some parameters are neither defined as set value nor random value
  all_parms <- c(parm_random,parm_constante)
  full_list_parms <- seq_len(length(input_prosail_names))
  set_to_default <- c()
  if (length(setdiff(full_list_parms,all_parms))>0){
    message('WARNING: elements are neither defined as parameters to sample')
    message('nor parameters to set: ')
    print(input_prosail_names[setdiff(full_list_parms,all_parms)])
    message('These parameters will be set to their default value')
    print(default[setdiff(full_list_parms,all_parms)])
    set_to_default <- setdiff(full_list_parms,all_parms)
  }

  # define input_prosail # 1 default value
  if (length(set_to_default)>0){
    for (i in set_to_default) input_prosail[[i]] <- default[,i]
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
    if (type_distrib[[sel]] == 'Uniform'){
      input_prosail[[sel]] <- array(runif(nb_samples,
                                          min = minval[1,i],
                                          max = maxval[1,i]),
                                    dim = c(nb_samples,1))
      # if Gaussian distribution
    } else if (type_distrib[[sel]] == 'Gaussian'){
      set.seed(42)
      input_prosail[[sel]] <- truncnorm::rtruncnorm(n = nb_samples,
                                                    a = minval[[sel]],
                                                    b = maxval[[sel]],
                                                    mean = mean[[sel]],
                                                    sd = sd[[sel]])
      input_prosail[[sel]] <- array(input_prosail[[sel]],dim = c(nb_samples,1))
    }
  }
  return(input_prosail)
}
