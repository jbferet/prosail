#' function identifying for which parameters prior information is known
#'
#' @param PriorInfoMean list. mean value expected for a list of parameters
#' @param PriorInfoSD list. dtandard deviation to the mean expected for a list
#' of parameters
#'
#' @return fc estimates of the parameters
#' @export
#'
which_parm_prior <- function(PriorInfoMean, PriorInfoSD) {

  Parms2Prior <- c()
  listParms <- c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'PROT', 'CBC', 'N',
                 'alpha', 'LIDFa', 'LIDFb', 'lai', 'q', 'tts', 'tto', 'psi',
                 'psoil')
  PriorInfoMean_Update <- PriorInfoSD_Update <- list()
  for (parm in listParms){
    if (parm %in% names(PriorInfoMean) & parm %in% names(PriorInfoSD)){
      Parms2Prior <- c(Parms2Prior,parm)
      PriorInfoMean_Update[[parm]] <- PriorInfoMean[[parm]]
      PriorInfoSD_Update[[parm]] <- PriorInfoSD[[parm]]
    }
  }
  PriorInfoMean_Update <- data.frame(PriorInfoMean_Update)
  PriorInfoSD_Update <- data.frame(PriorInfoSD_Update)
  return(list('Parms2Prior' = Parms2Prior,
              'PriorInfoMean' = PriorInfoMean_Update,
              'PriorInfoSD' = PriorInfoSD_Update))
}
