#' Value of the cost criterion to minimize during PROSAIL inversion
#'
#' @param refl_mes numeric. Measured surface reflectance
#' @param refl_mod numeric. Simulated surface reflectance
#' @param xprior list. values of the parameters for which prior information is
#'  provided
#' @param prior_info list. prior mean, sd and weight of parameters defined as
#' xprior modulate its importance
#'
#' @return cost
#' @export
cost_function_rmse_prosail  <- function(refl_mes, refl_mod, xprior,
                                        prior_info = NULL){

  fc <- sqrt(sum((refl_mes-refl_mod)**2)/length(refl_mes))
  if (!is.null(prior_info))
    fc <- fc +
      prior_info$weight_prior*mean(as.numeric((xprior-prior_info$mean)/prior_info$sd)**2)
  return(fc)
}


#' @rdname prosail-deprecated
#' @export
CostVal_RMSE_PROSAIL  <- function(brfMES, brfMOD, xprior, PriorInfoMean = NULL,
                                  PriorInfoSD = NULL, WeightPrior = 0.01){
  .Deprecated("cost_function_rmse_prosail")
  prior_info <- list('mean' = PriorInfoMean,
                     'SD' = PriorInfoSD,
                     'weight_prior' = WeightPrior)
  cost_function_rmse_prosail(brfMES, brfMOD, xprior, prior_info)
}
