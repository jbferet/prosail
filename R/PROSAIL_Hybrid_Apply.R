#' This function applies the regression models trained with prosail_hybrid_train
#'
#' @param regression_models list. regression models produced by
#' prosail_hybrid_train
#' @param refl numeric. LUT of BRF used for training
#' @param progressBar boolean. should progressbar be displayed?
#'
#' @return HybridRes list. Estimated values corresponding to refl. Includes
#' - MeanEstimate = mean value for the ensemble regression model
#' - StdEstimate = std value for the ensemble regression model
#' @importFrom stats predict
#' @importFrom matrixStats rowSds
#' @importFrom progress progress_bar
#' @importFrom methods is
#' @export

prosail_hybrid_apply <- function(regression_models, refl, progressBar = FALSE){

  # make sure refl is right dimensions
  refl <- t(refl)
  if (inherits(regression_models[[1]], what = 'liquidSVM')){
    nbFeatures <- regression_models[[1]]$dim
  } else {
    nbFeatures <- ncol(regression_models[[1]]$trainingData) - 1
  }
  if (!ncol(refl)==nbFeatures & nrow(refl)==nbFeatures){
    refl <- t(refl)
  }
  nb_bagg <- length( regression_models)
  EstimatedVal <- list()
  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Applying SVR models [:bar] :percent in :elapsed",
      total = nb_bagg, clear = FALSE, width= 100)
  }
  for (i in seq_len(nb_bagg)){
    EstimatedVal[[i]] <- predict(regression_models[[i]], refl)
    if (progressBar == TRUE) pb$tick()
  }
  EstimatedVal <- do.call(cbind,EstimatedVal)
  MeanEstimate <- rowMeans(EstimatedVal)
  StdEstimate <- rowSds(EstimatedVal)
  return(list("MeanEstimate" = MeanEstimate,
              "StdEstimate" = StdEstimate))
}
