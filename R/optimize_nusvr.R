#' This function optimizes kernlab nu svr model using bayesian optimization
#'
#' @param lut numeric. LUT of BRF used for training
#' @param options numeric. biophysical parameter corresponding to refl
#'
#' @return list of optimized hyperparameters
#' @importFrom mlr makeRegrTask makeTuneControlRandom makeLearner makeResampleDesc tuneParams
#' @importFrom ParamHelpers makeParamSet makeNumericParam
#' @export

optimize_nusvr <- function(lut, options = NULL){

  options <- set_options_prosail(fun = 'optimize_nusvr',
                                 options = options)
  task <- mlr::makeRegrTask(id = deparse(substitute(lut)), data = lut,
                            target = "target")
  ps <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("C", lower = options$C$lower,
                                   upper = options$C$upper,
                                   trafo = function(x) log(x)),
    ParamHelpers::makeNumericParam("nu", lower = options$nu$lower,
                                   upper = options$nu$upper),
    ParamHelpers::makeNumericParam("sigma", lower = options$sigma$lower,
                                   upper = options$sigma$upper,
                                   trafo = function(x) log(x)))

  ctrl <- mlr::makeTuneControlRandom(maxit = options$maxit)
  learner <- mlr::makeLearner(cl = "regr.ksvm", type = "nu-svr")
  rdesc <- mlr::makeResampleDesc("CV", iters = 5)
  tune <- mlr::tuneParams(learner = learner, task = task, par.set = ps,
                          control = ctrl, show.info = F, resampling = rdesc)
  return(tune$x)
}
