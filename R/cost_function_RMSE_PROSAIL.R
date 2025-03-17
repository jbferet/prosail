#' Value of the cost criterion to minimize during PROSAIL inversion
#' @param brfMES numeric. Measured bidirectional reflectance
#' @param brfMOD numeric. Simulated bidirectional reflectance
#' @param xprior list. values of the parameters for which prior information is
#'  provided
#' @param PriorInfoMean list. prior mean value of parameters defined as xprior
#' @param PriorInfoSD list. prior standard dev of parameters defined as xprior
#' @param WeightPrior numeric. Weight to be applied on prior information to
#' modulate its importance
#'
#' @return res list. Includes Parms2Estimate, Parms2Set,
#' InitialGuess, LowerBound, UpperBound, ParmSet, InVar
#' @export
cost_function_RMSE_PROSAIL  <- function(brfMES, brfMOD, xprior,
                                        PriorInfoMean = NULL,
                                        PriorInfoSD = NULL,
                                        WeightPrior = 0.01){

  fc <- sqrt(sum((brfMES-brfMOD)**2)/length(brfMES))
  if (!is.null(PriorInfoMean))
    fc <- fc +
      WeightPrior*mean(as.numeric((xprior-PriorInfoMean)/PriorInfoSD)**2)
  return(fc)
}
