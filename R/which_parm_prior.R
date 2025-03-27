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

  Parms2Prior <- c()
  listParms <- c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'PROT', 'CBC', 'N',
                 'alpha', 'LIDFa', 'LIDFb', 'lai', 'q', 'tts', 'tto', 'psi',
                 'psoil')
  prior_mean_update <- prior_sd_update <- list()
  for (parm in listParms){
    if (parm %in% names(prior_mean) & parm %in% names(prior_sd)){
      Parms2Prior <- c(Parms2Prior,parm)
      prior_mean_update[[parm]] <- prior_mean[[parm]]
      prior_sd_update[[parm]] <- prior_sd[[parm]]
    }
  }
  prior_mean_update <- data.frame(prior_mean_update)
  prior_sd_update <- data.frame(prior_sd_update)
  return(list('Parms2Prior' = Parms2Prior,
              'mean' = prior_mean_update,
              'sd' = prior_sd_update))
}
