#' function identifying which parameters should be estimated during inversion
#'
#' @param InitialGuess list. Initial values during inversion
#' @param LowerBound list. Lower bound values during inversion
#' @param UpperBound list. Upper bound values during inversion
#' @param ParmSet list. Values for variables out of inversion
#'
#' @return res list. includes
#' - Parms2Estimate
#' - Parms2Set
#' - InitialGuess
#' - LowerBound
#' - UpperBound
#' - ParmSet
#' - InVar
#' fc estimates of the parameters
#' @export
#'
which_parms_to_invert <- function(InitialGuess, LowerBound,
                                   UpperBound, ParmSet) {

  # define all parameters which can be assessed through iterativbe optimization
  InVar <- data.frame('CHL' = 0,'CAR' = 0,'ANT' = 0,'BROWN' = 0,'EWT' = 0,
                      'LMA' = 0,'PROT' = 0,'CBC' = 0,'N' = 0,'alpha' = 40,
                      'LIDFa' = 0,'LIDFb' = 0,'lai' = 0,'q' = 0,
                      'tts' = 0,'tto' = 0,'psi' = 0,'psoil' = 0)
  AllParms <- names(InVar)
  Parms2Estimate <- ParmSet_Final <- c()
  InitialGuess_Update <- LowerBound_Update <- UpperBound_Update <-
    ParmSet_Update <- list()

  # set parameters to user value defined in ParmSet
  for (parm in AllParms){
    if (parm %in% names(InitialGuess) & !parm %in% names(ParmSet)){
      Parms2Estimate <- c(Parms2Estimate,parm)
      InitialGuess_Update[[parm]] <- InitialGuess[[parm]]
      LowerBound_Update[[parm]] <- LowerBound[[parm]]
      UpperBound_Update[[parm]] <- UpperBound[[parm]]
    } else if (!parm %in% names(InitialGuess) & parm %in% names(ParmSet)){
      ParmSet_Final <- c(ParmSet_Final,parm)
      ParmSet_Update[[parm]] <- ParmSet[[parm]]
    }
  }
  InitialGuess_Update <- data.frame(InitialGuess_Update)
  LowerBound_Update <- data.frame(LowerBound_Update)
  UpperBound_Update <- data.frame(UpperBound_Update)
  ParmSet_Update <- data.frame(ParmSet_Update)
  return(list('Parms2Estimate' = Parms2Estimate,
              'Parms2Set' = ParmSet_Final,
              'InitialGuess' = InitialGuess_Update,
              'LowerBound' = LowerBound_Update,
              'UpperBound' = UpperBound_Update,
              'ParmSet' = ParmSet_Update,
              'InVar' = InVar))
}
