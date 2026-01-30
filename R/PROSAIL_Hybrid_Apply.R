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
#' @importFrom kernlab xmatrix
#' @export

prosail_hybrid_apply <- function(regression_models, refl, progressBar = FALSE){

  # make sure refl is right dimensions
  refl <- t(refl)
  if (inherits(regression_models[[1]], what = 'liquidSVM')){
    nb_features <- regression_models[[1]]$dim
  } else if (inherits(regression_models[[1]], what = 'train')){
    nb_features <- ncol(regression_models[[1]]$trainingData) - 1
  } else if (inherits(regression_models[[1]], what = 'ksvm')){
    nb_features <- ncol(kernlab::xmatrix(regression_models[[1]]))
  }
  if (!ncol(refl)==nb_features & nrow(refl)==nb_features)
    refl <- t(refl)
  nb_bagg <- length( regression_models)
  estimated_val <- list()
  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Applying SVR models [:bar] :percent in :elapsed",
      total = nb_bagg, clear = FALSE, width= 100)
  }
  for (i in seq_len(nb_bagg)){
    if (!inherits(regression_models[[i]], what = 'ksvm')){
      estimated_val[[i]] <- stats::predict(object = regression_models[[i]],
                                           refl)
    } else if (inherits(regression_models[[i]], what = 'ksvm')){
      estimated_val[[i]] <- kernlab::predict(object = regression_models[[i]],
                                             refl)
    }
    if (progressBar == TRUE)
      pb$tick()
  }
  estimated_val <- do.call(cbind,estimated_val)
  mean_estimate <- rowMeans(estimated_val)
  sd_estimate <- rowSds(estimated_val)
  return(list("MeanEstimate" = mean_estimate,
              "StdEstimate" = sd_estimate))
}


#' @rdname prosail-deprecated
#' @export
PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl, progressBar = FALSE){
  .Deprecated("prosail_hybrid_apply")
  prosail_hybrid_apply(RegressionModels,Refl, progressBar)
}

