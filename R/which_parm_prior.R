#' function identifying for which parameters prior information is known
#'
#' @param prior_mean list. mean value expected for a list of parameters
#' @param prior_sd list. standard deviation expected for a list
#' of parameters
#'
#' @return fc estimates of the parameters
#' @export
#'
which_parm_prior <- function(prior_mean, prior_sd) {

  parms_to_prior <- c()
  list_parms <- c('chl', 'car', 'ant', 'brown', 'ewt', 'lma', 'prot', 'cbc',
                  'n_struct', 'alpha', 'lidf_a', 'lidf_b', 'lai', 'hotspot', 'tts',
                  'tto', 'psi', 'psoil')
  prior_mean_update <- prior_sd_update <- list()
  for (parm in list_parms){
    if (parm %in% names(prior_mean) & parm %in% names(prior_sd)){
      parms_to_prior <- c(parms_to_prior,parm)
      prior_mean_update[[parm]] <- prior_mean[[parm]]
      prior_sd_update[[parm]] <- prior_sd[[parm]]
    }
  }
  prior_mean_update <- data.frame(prior_mean_update)
  prior_sd_update <- data.frame(prior_sd_update)
  return(list('parms_to_prior' = parms_to_prior,
              'mean' = prior_mean_update,
              'sd' = prior_sd_update))
}
