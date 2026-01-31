#' This function trains a regression model to estimate a set of variables
#' from spectral data
#'
#' @param refl_lut numeric. LUT of surface reflectance used for training
#' @param input_variables numeric. biophysical parameter corresponding to refl
#' @param nb_bagg numeric. nb of individual subsets generated from refl_lut
#' @param replacement Boolean. subsets generated with / without replacement
#' @param method character. which machine learning regression method used?
#' default = SVM with liquidSVM.
#' set nu-svr for RBF SVM with kernlab package
#' set svmRadial or svmLinear from SVM with caret package
#' More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#' @param progressBar boolean. should progressbar be displayed?
#' @param options list
#'
#' @return models_mlr list. ML regression models trained for the retrieval of
#' input_variables based on refl_lut
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom simsalapar tryCatch.W.E
#' @importFrom caret train trainControl
#' @importFrom kernlab ksvm
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @export

prosail_hybrid_train <- function(refl_lut, input_variables, nb_bagg = 20,
                                 replacement = FALSE, method = 'liquidSVM',
                                 verbose = FALSE, progressBar = FALSE,
                                 options = NULL){
  is_liquidSVM_available <- system.file(package = "liquidSVM")
  if (method == 'liquidSVM' & is_liquidSVM_available=='')
    method <- 'nu-svr'
  x <- y <- ymean <- ystdmin <- ystdmax <- NULL
  # split the LUT into nb_bagg subsets
  nb_samples <- length(input_variables)
  if (dim(refl_lut)[2]==nb_samples)
    refl_lut <- t(refl_lut)
  # if subsets are generated from refl_lut with replacement
  if (replacement==TRUE){
    subsets <- list()
    samples_per_run <- round(nb_samples/nb_bagg)
    for (run in seq_len(nb_bagg))
      subsets[[run]] <- sample(seq_len(nb_samples), samples_per_run,
                               replace = TRUE)
    # if subsets are generated from refl_lut without replacement
  } else if (replacement==FALSE){
    subsets <- split(sample(seq_len(nb_samples)),seq_len(nb_bagg))
  }

  # run training for each subset
  models_mlr <- list()

  if (progressBar == TRUE)
    pb <- progress_bar$new(
      format = "Training SVR [:bar] :percent in :elapsedfull , eta = :eta",
      total = nb_bagg, clear = FALSE, width= 100)

  refl_lut <- data.frame(refl_lut)
  for (i in seq_len(nb_bagg)){
    training_set <- list()
    training_set$X <- refl_lut %>% dplyr::slice(subsets[i][[1]])
    training_set$Y <- input_variables[subsets[i][[1]]]

    # if using caret
    control <- caret::trainControl(method="cv", number = 5)
    if (is.null(colnames(training_set$X)))
      colnames(training_set$X) <- paste('var', seq_len(ncol(training_set$X)),
                                        sep = '_')
    target <- matrix(training_set$Y,ncol = 1)
    if (is.null(colnames(target)))
      colnames(target) <- 'target'
    training_data <- cbind(target,training_set$X)

    if (method == 'liquidSVM'){
      r1 <- tryCatch.W.E(tuned_model <- liquidSVM::svmRegression(training_set$X,
                                                                 training_set$Y))
      if (!is.null(r1$warning)){
        msg <- r1$warning$message
        val_gamma <- stringr::str_split(string = msg,
                                        pattern = 'gamma=')[[1]][2]
        val_lambda <- stringr::str_split(string = msg,
                                         pattern = 'lambda=')[[1]][2]
        if (!is.na(as.numeric(val_gamma))){
          if (verbose==TRUE)
          { message('Adjusting gamma accordingly')}
          val_gamma <- as.numeric(val_gamma)
          tuned_model <- liquidSVM::svmRegression(training_set$X,
                                                  training_set$Y,
                                                  min_gamma = val_gamma)
        }
        if (!is.na(as.numeric(val_lambda))){
          if (verbose==TRUE) { message('Adjusting lambda accordingly')}
          val_lambda <- as.numeric(val_lambda)
          tuned_model <- liquidSVM::svmRegression(training_set$X,
                                                  training_set$Y,
                                                  min_lambda = val_lambda)
        }
      }
    } else if (method %in% c('nu-svr')){
      # Define parameter grid
      C_values <- c(0.1, 1, 10, 100)
      nu_values <- c(0.01, 0.1, 0.25, 0.5, 0.75)
      sigma_values <- c(0.001, 0.01, 0.1, 1, 10)
      best_model <- NULL
      min_err <- Inf
      lut <- data.frame(training_set$X, 'target' = training_set$Y)
      for (C_val in C_values) {
        for (nu_val in nu_values) {
          for (sigma_val in sigma_values) {
            model <- kernlab::ksvm(
              target ~ .,
              data = lut,
              type = "nu-svr",
              kernel = "rbfdot",
              kpar = list(sigma = sigma_val),
              C = C_val,
              nu = nu_val,
              cross = 5  # Cross-validation folds
            )
            err_cv <- model@cross  # Mean squared error from CV
            if (err_cv < min_err) {  # Negative because it's an error
              min_err <- err_cv
              best_model <- model
              optim_parms <- list('C' = C_val, 'nu' =nu_val, 'sigma' = sigma_val)
            }
          }
        }
      }
      tuned_model <- kernlab::ksvm(
        target ~ ., data = lut,
        type="nu-svr", kernel="rbfdot",
        kpar <- list(optim_parms$sigma),
        C <- optim_parms$C,
        nu <- optim_parms$nu)

      # # optimize hyperparameters based on bayesian optimization
      # lut <- data.frame(training_set$X, 'target' = training_set$Y)
      # optim_parms <- optimize_nusvr(lut = lut, options = options)
      # # adjust model with tuned parms
      # tuned_model <- kernlab::ksvm(
      #   target ~ ., data = lut,
      #   type="nu-svr", kernel="rbfdot",
      #   kpar <- list(optim_parms$sigma),
      #   C <- optim_parms$C,
      #   nu <- optim_parms$nu)

    } else if (method %in% c('svmRadial', 'svmLinear')){
      if (method =='svmRadial')
        tuneGrid <- expand.grid(C = exp(seq(-10, 0)),
                                sigma = exp(seq(-10, 0)))
      if (method =='svmLinear')
        tuneGrid <- expand.grid(C = exp(seq(-10, 0)))
      tuned_model <- caret::train(target ~ .,
                                  data = training_data,
                                  method = method,
                                  preProcess = c("center", "scale"),
                                  trControl = control,
                                  tuneGrid = tuneGrid)
    } else if (method == 'gaussprLinear'){
      tuned_model <- caret::train(target ~ .,
                                  data = training_data,
                                  method = method,
                                  preProcess = c("center", "scale"),
                                  trControl = control)
    } else if (method == 'rf'){
      seed <- 7
      set.seed(seed)
      mtry <- ncol(training_set$X)/3
      tuned_model <- train(target~.,
                           data = training_data,
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
      #   tuned_model <- caret::train(
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
      #   tuned_model <- train(target~.,
      #                       data = training_data,
      #                       method = method,
      #                       metric = 'RMSE',
      #                       preProcess = c("center", "scale"),
      #                       tuneGrid = grid,
      #                       gamma = 0.5,
      #                       trControl = control)
    }
    models_mlr[[i]] <- tuned_model
    if (progressBar == TRUE) pb$tick()
  }
  return(models_mlr)
}


#' @rdname prosail-deprecated
#' @export

PROSAIL_Hybrid_Train <- function(BRF_LUT, InputVar, nbEnsemble = 20,
                                 WithReplacement = FALSE,
                                 method = 'liquidSVM',
                                 verbose = FALSE, progressBar = FALSE){
  .Deprecated("prosail_hybrid_train")
  prosail_hybrid_train(BRF_LUT, InputVar, nbEnsemble, WithReplacement,
                       method, verbose, progressBar)
}

