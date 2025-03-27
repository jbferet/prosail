#' Value of the cost criterion to minimize during PROSAIL inversion
#'
#' @param brf_mes numeric. Measured bidirectional reflectance
#' @param brf_mod numeric. Simulated bidirectional reflectance
#' @param xprior list. values of the parameters for which prior information is
#'  provided
#' @param prior_info list. prior mean, sd and weight of parameters defined as
#' xprior modulate its importance
#'
#' @return cost
#' @export
cost_function_RMSE_PROSAIL  <- function(brf_mes, brf_mod, xprior,
                                        prior_info = NULL){

  fc <- sqrt(sum((brf_mes-brf_mod)**2)/length(brf_mes))
  if (!is.null(prior_info))
    fc <- fc +
      prior_info$WeightPrior*mean(as.numeric((xprior-prior_info$mean)/prior_info$sd)**2)
  return(fc)
}
