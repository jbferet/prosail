#' function identifying which parameters should be estimated during inversion
#'
#' @param initialization list. Initial values during inversion
#' @param lower_bound list. Lower bound values during inversion
#' @param upper_bound list. Upper bound values during inversion
#' @param parm_set list. Values for variables out of inversion
#'
#' @return res list. includes
#' - Parms2Estimate
#' - Parms2Set
#' - initialization
#' - lower_bound
#' - upper_bound
#' - parm_set
#' - InVar
#' fc estimates of the parameters
#' @export
#'
which_parms_to_invert <- function(initialization, lower_bound,
                                  upper_bound, parm_set) {

  # define all parameters which can be assessed through iterativbe optimization
  input_prosail <- data.frame('CHL' = 0, 'CAR' = 0, 'ANT' = 0, 'BROWN' = 0,
                              'EWT' = 0, 'LMA' = 0, 'PROT' = 0, 'CBC' = 0,
                              'N' = 0, 'alpha' = 40, 'LIDFa' = 0, 'LIDFb' = 0,
                              'lai' = 0, 'q' = 0, 'tts' = 0, 'tto' = 0,
                              'psi' = 0, 'psoil' = 0)
  AllParms <- names(input_prosail)
  Parms2Estimate <- parm_set_Final <- c()
  initialization_Update <- lower_bound_Update <- upper_bound_Update <-
    parm_set_Update <- list()

  # set parameters to user value defined in parm_set
  for (parm in AllParms){
    if (parm %in% names(initialization) & !parm %in% names(parm_set)){
      Parms2Estimate <- c(Parms2Estimate,parm)
      initialization_Update[[parm]] <- initialization[[parm]]
      lower_bound_Update[[parm]] <- lower_bound[[parm]]
      upper_bound_Update[[parm]] <- upper_bound[[parm]]
    } else if (!parm %in% names(initialization) & parm %in% names(parm_set)){
      parm_set_Final <- c(parm_set_Final,parm)
      parm_set_Update[[parm]] <- parm_set[[parm]]
    }
  }
  initialization_Update <- data.frame(initialization_Update)
  lower_bound_Update <- data.frame(lower_bound_Update)
  upper_bound_Update <- data.frame(upper_bound_Update)
  parm_set_Update <- data.frame(parm_set_Update)
  return(list('Parms2Estimate' = Parms2Estimate,
              'Parms2Set' = parm_set_Final,
              'initialization' = initialization_Update,
              'lower_bound' = lower_bound_Update,
              'upper_bound' = upper_bound_Update,
              'parm_set' = parm_set_Update,
              'InVar' = input_prosail))
}
