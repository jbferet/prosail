#' function identifying which parameters should be estimated during inversion
#'
#' @param initialization list. Initial values during inversion
#' @param lower_bound list. Lower bound values during inversion
#' @param upper_bound list. Upper bound values during inversion
#' @param parm_set list. Values for variables out of inversion
#'
#' @return res list. includes
#' - parms_to_estimate
#' - parms_to_set
#' - initialization
#' - lower_bound
#' - upper_bound
#' - parm_set
#' - input_prosail
#' fc estimates of the parameters
#' @export
#'
which_parms_to_invert <- function(initialization, lower_bound,
                                  upper_bound, parm_set) {

  # define all parameters which can be assessed through iterativbe optimization
  input_prosail <- data.frame('chl' = 0, 'car' = 0, 'ant' = 0, 'brown' = 0,
                              'ewt' = 0, 'lma' = 0, 'prot' = 0, 'cbc' = 0,
                              'n_struct' = 0, 'alpha' = 40, 'lidf_a' = 0,
                              'lidf_b' = 0, 'lai' = 0, 'hotspot' = 0, 'tts' = 0,
                              'tto' = 0, 'psi' = 0, 'psoil' = 0)
  all_parms <- names(input_prosail)
  parms_to_estimate <- parm_set_final <- c()
  initialization_update <- lower_bound_update <- upper_bound_update <-
    parm_set_update <- list()

  # set parameters to user value defined in parm_set
  for (parm in all_parms){
    if (parm %in% names(initialization) & !parm %in% names(parm_set)){
      parms_to_estimate <- c(parms_to_estimate,parm)
      initialization_update[[parm]] <- initialization[[parm]]
      lower_bound_update[[parm]] <- lower_bound[[parm]]
      upper_bound_update[[parm]] <- upper_bound[[parm]]
    } else if (!parm %in% names(initialization) & parm %in% names(parm_set)){
      parm_set_final <- c(parm_set_final,parm)
      parm_set_update[[parm]] <- parm_set[[parm]]
    }
  }
  initialization_update <- data.frame(initialization_update)
  lower_bound_update <- data.frame(lower_bound_update)
  upper_bound_update <- data.frame(upper_bound_update)
  parm_set_update <- data.frame(parm_set_update)
  return(list('parms_to_estimate' = parms_to_estimate,
              'parms_to_set' = parm_set_final,
              'initialization' = initialization_update,
              'lower_bound' = lower_bound_update,
              'upper_bound' = upper_bound_update,
              'parm_set' = parm_set_update,
              'input_prosail' = input_prosail))
}
