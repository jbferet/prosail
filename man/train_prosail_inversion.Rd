% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_PROSAIL_HybridInversion.R
\name{train_prosail_inversion}
\alias{train_prosail_inversion}
\title{This function performs full training for hybrid inversion using SVR with
values for default parameters}
\usage{
train_prosail_inversion(
  InputPROSAIL = NULL,
  BRF_LUT = NULL,
  atbd = FALSE,
  GeomAcq = NULL,
  Codist_LAI = TRUE,
  minval = NULL,
  maxval = NULL,
  TypeDistrib = NULL,
  GaussianDistrib = NULL,
  ParmSet = NULL,
  SAILversion = "4SAIL",
  BrownLOP = NULL,
  nbSamples = 2000,
  nbSamplesPerRun = 100,
  nbModels = 20,
  Replacement = TRUE,
  Parms2Estimate = "lai",
  Bands2Select = NULL,
  NoiseLevel = NULL,
  SRF = NULL,
  SpecPROSPECT = NULL,
  SpecSOIL = NULL,
  SpecATM = NULL,
  Path_Results = "./",
  FigPlot = FALSE,
  method = "liquidSVM",
  verbose = FALSE
)
}
\arguments{
\item{InputPROSAIL}{list. user-defined list of input parameters to be used to produce a training LUT}

\item{BRF_LUT}{list. user-defined BRF LUT used to run the hybrid inversion}

\item{atbd}{boolean. should input parameter distribution from ATBD be applied ?}

\item{GeomAcq}{list. geometry of acquisiton. list should contain min and max values for tts, tto and psi}

\item{Codist_LAI}{boolean. set TYRUE if codistribution with LAI accounted for}

\item{minval}{list. minimum value for input parameters sampled to produce a training LUT}

\item{maxval}{list. maximum value for input parameters sampled to produce a training LUT}

\item{TypeDistrib}{list. Type of distribution. Either 'Uniform' or 'Gaussian'}

\item{GaussianDistrib}{list. Mean value and STD corresponding to the parameters sampled with gaussian distribution}

\item{ParmSet}{list. list of input parameters set to a specific value}

\item{SAILversion}{character. Either 4SAIL or 4SAIL2}

\item{BrownLOP}{character. Either 4SAIL or 4SAIL2}

\item{nbSamples}{numeric. number of samples in training LUT}

\item{nbSamplesPerRun}{numeric. number of training sample per individual regression model}

\item{nbModels}{numeric. number of individual models to be run for ensemble}

\item{Replacement}{bolean. is there replacement in subsampling?}

\item{Parms2Estimate}{list. list of input parameters to be estimated}

\item{Bands2Select}{list. list of bands used for regression for each input parameter}

\item{NoiseLevel}{list. list of noise value added to reflectance (defined per input parm)}

\item{SRF}{list. Spectral response function}

\item{SpecPROSPECT}{list. Includes optical constants required for PROSPECT}

\item{SpecSOIL}{list. Includes either dry soil and wet soil, or a unique soil sample if the psoil parameter is not inverted}

\item{SpecATM}{list. Includes direct and diffuse radiation for clear conditions}

\item{Path_Results}{character. path for results}

\item{FigPlot}{boolean. Set TRUE to get scatterplot of estimated biophysical variable during training step}

\item{method}{character. which machine learning regression method should be used?
default = SVM with liquidSVM. svmRadial and svmLinear from caret package also implemented. More to come}

\item{verbose}{boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training}
}
\value{
modelsSVR list. regression models trained for the retrieval of InputVar based on BRF_LUT
}
\description{
This function performs full training for hybrid inversion using SVR with
values for default parameters
}
