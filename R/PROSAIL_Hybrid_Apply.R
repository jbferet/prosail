#' This function applies the regression models trained with PROSAIL_Hybrid_Train
#'
#' @param RegressionModels list. regression models produced by
#' PROSAIL_Hybrid_Train
#' @param Refl numeric. LUT of BRF used for training
#' @param progressBar boolean. should progressbar be displayed?
#'
#' @return HybridRes list. Estimated values corresponding to Refl. Includes
#' - MeanEstimate = mean value for the ensemble regression model
#' - StdEstimate = std value for the ensemble regression model
#' @importFrom stats predict
#' @importFrom matrixStats rowSds
#' @importFrom progress progress_bar
#' @importFrom methods is
#' @export

PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl, progressBar = FALSE){

  # make sure Refl is right dimensions
  Refl <- t(Refl)
  if (inherits(RegressionModels[[1]], what = 'liquidSVM')){
    nbFeatures <- RegressionModels[[1]]$dim
  } else {
    nbFeatures <- ncol(RegressionModels[[1]]$trainingData) - 1
  }
  if (!ncol(Refl)==nbFeatures & nrow(Refl)==nbFeatures){
    Refl <- t(Refl)
  }
  nb_bagg <- length( RegressionModels)
  EstimatedVal <- list()
  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Applying SVR models [:bar] :percent in :elapsed",
      total = nb_bagg, clear = FALSE, width= 100)
  }
  for (i in seq_len(nb_bagg)){
    EstimatedVal[[i]] <- predict(RegressionModels[[i]], Refl)
    if (progressBar == TRUE) pb$tick()
  }
  EstimatedVal <- do.call(cbind,EstimatedVal)
  MeanEstimate <- rowMeans(EstimatedVal)
  StdEstimate <- rowSds(EstimatedVal)
  return(list("MeanEstimate" = MeanEstimate,
              "StdEstimate" = StdEstimate))
}
