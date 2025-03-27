#' This function trains a support vector regression for a set of variables based on spectral data
#'
#' @param BRF_LUT numeric. LUT of BRF used for training
#' @param input_variables numeric. biophysical parameter corresponding to reflectance
#' @param nb_bagg numeric. nb of individual subsets generated from BRF_LUT
#' @param replacement Boolean. subsets generated with / without replacement
#' @param method character. which machine learning regression method used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package
#' also implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#' @param progressBar boolean. should progressbar be displayed?
#'
#' @return modelsMLR list. ML regression models trained for the retrieval of
#' input_variables based on BRF_LUT
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom simsalapar tryCatch.W.E
#' @importFrom caret train trainControl
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @export

PROSAIL_Hybrid_Train <- function(BRF_LUT, input_variables, nb_bagg = 20,
                                 replacement = FALSE, method = 'liquidSVM',
                                 verbose = FALSE, progressBar = FALSE){

  x <- y <- ymean <- ystdmin <- ystdmax <- NULL
  # split the LUT into nb_bagg subsets
  nb_samples <- length(input_variables)
  if (dim(BRF_LUT)[2]==nb_samples)
    BRF_LUT <- t(BRF_LUT)
  # if subsets are generated from BRF_LUT with replacement
  if (replacement==TRUE){
    Subsets <- list()
    samples_per_run <- round(nb_samples/nb_bagg)
    for (run in seq_len(nb_bagg))
      Subsets[[run]] <- sample(seq_len(nb_samples), samples_per_run,
                               replace = TRUE)
    # if subsets are generated from BRF_LUT without replacement
  } else if (replacement==FALSE){
    Subsets <- split(sample(seq_len(nb_samples)),seq_len(nb_bagg))
  }

  # run training for each subset
  modelsMLR <- predictedYAll <- tunedModelYAll <- list()

  if (progressBar == TRUE){
    pb <- progress_bar$new(
      format = "Training SVR [:bar] :percent in :elapsedfull , eta = :eta",
      total = nb_bagg, clear = FALSE, width= 100)
  }
  BRF_LUT <- data.frame(BRF_LUT)
  for (i in seq_len(nb_bagg)){
    TrainingSet <- list()
    TrainingSet$X <- BRF_LUT %>% dplyr::slice(Subsets[i][[1]])
    TrainingSet$Y <- input_variables[Subsets[i][[1]]]

    # if using caret
    control <- caret::trainControl(method="cv", number = 5)
    if (is.null(colnames(TrainingSet$X)))
      colnames(TrainingSet$X) <- paste('var', seq_len(ncol(TrainingSet$X)),
                                       sep = '_')
    target <- matrix(TrainingSet$Y,ncol = 1)
    if (is.null(colnames(target)))
      colnames(target) <- 'target'
    TrainingData <- cbind(target,TrainingSet$X)

    if (method == 'liquidSVM'){
      # liquidSVM
      r1 <- tryCatch.W.E(tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                                TrainingSet$Y))
      # tunedModel <- liquidSVM::svmRegression(TrainingSet$X, TrainingSet$Y)
      if (!is.null(r1$warning)){
        Msg <- r1$warning$message
        ValGamma <- stringr::str_split(string = Msg,
                                       pattern = 'gamma=')[[1]][2]
        ValLambda <- stringr::str_split(string = Msg,
                                        pattern = 'lambda=')[[1]][2]
        if (!is.na(as.numeric(ValGamma))){
          if (verbose==TRUE)
          { message('Adjusting Gamma accordingly')}
          ValGamma <- as.numeric(ValGamma)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                 TrainingSet$Y,
                                                 min_gamma = ValGamma)
        }
        if (!is.na(as.numeric(ValLambda))){
          if (verbose==TRUE) { message('Adjusting Lambda accordingly')}
          ValLambda <- as.numeric(ValLambda)
          tunedModel <- liquidSVM::svmRegression(TrainingSet$X,
                                                 TrainingSet$Y,
                                                 min_lambda = ValLambda)
        }
      }
    } else if (method %in% c('svmRadial', 'svmLinear')){
      if (method =='svmRadial')
        tuneGrid <- expand.grid(C = exp(seq(-10, 0)),
                                sigma = exp(seq(-10, 0)))
      if (method =='svmLinear')
        tuneGrid <- expand.grid(C = exp(seq(-10, 0)))
      tunedModel <- caret::train(target ~ .,
                                 data = TrainingData,
                                 method = method,
                                 preProcess = c("center", "scale"),
                                 trControl = control,
                                 tuneGrid = tuneGrid)
    } else if (method == 'gaussprLinear'){
      tunedModel <- caret::train(target ~ .,
                                 data = TrainingData,
                                 method = method,
                                 preProcess = c("center", "scale"),
                                 trControl = control)
    } else if (method == 'rf'){
      seed <- 7
      set.seed(seed)
      mtry <- ncol(TrainingSet$X)/3
      tunedModel <- train(target~.,
                          data = TrainingData,
                          method = method,
                          metric = 'RMSE',
                          preProcess = c("center", "scale"),
                          tuneLength = 15,
                          trControl = control)
      # } else if (method == 'xgbLinear'){
      #   grid_default <- expand.grid(nrounds = 100,
      #                               max_depth = 6,
      #                               eta = 0.3,
      #                               gamma = 0,
      #                               colsample_bytree = 1,
      #                               min_child_weight = 1,
      #                               subsample = 1)
      #   train_control <- caret::trainControl(method = "none",
      #                                        verboseIter = FALSE,
      #                                        allowParallel = TRUE )
      #
      #   tunedModel <- caret::train(
      #     x = trainMNX,
      #     y = trainMNY,
      #     trControl = train_control,
      #     tuneGrid = grid_default,
      #     method = "xgbTree",
      #     verbose = TRUE,
      #     nthreads = 4
      #   )
      #
      #
      #   grid <- expand.grid(nrounds = c(10,20), lambda = c(0.1),
      #                       alpha = c(1), eta = c(0.1))
      #   tunedModel <- train(target~.,
      #                       data = TrainingData,
      #                       method = method,
      #                       metric = 'RMSE',
      #                       preProcess = c("center", "scale"),
      #                       tuneGrid = grid,
      #                       gamma = 0.5,
      #                       trControl = control)
    }
    modelsMLR[[i]] <- tunedModel
    if (progressBar == TRUE) pb$tick()
  }
  return(modelsMLR)
}
