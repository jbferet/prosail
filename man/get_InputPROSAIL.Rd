% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_PROSAIL_HybridInversion.R
\name{get_InputPROSAIL}
\alias{get_InputPROSAIL}
\title{generate InputPROSAIL, following
\itemize{
\item either distribution defined in ATBD
\item or distribution defined in user
}}
\usage{
get_InputPROSAIL(
  atbd = FALSE,
  GeomAcq = NULL,
  Codist_LAI = TRUE,
  minval = NULL,
  maxval = NULL,
  TypeDistrib = NULL,
  GaussianDistrib = NULL,
  ParmSet = NULL,
  nbSamples = 2000,
  verbose = FALSE
)
}
\arguments{
\item{atbd}{boolean. should input parameter distribution from ATBD be applied ?}

\item{GeomAcq}{list. geometry of acquisiton. list should contain min and max values for tts, tto and psi}

\item{Codist_LAI}{boolean. set TYRUE if codistribution with LAI accounted for}

\item{minval}{list. minimum value for input parameters sampled to produce a training LUT}

\item{maxval}{list. maximum value for input parameters sampled to produce a training LUT}

\item{TypeDistrib}{list. Type of distribution. Either 'Uniform' or 'Gaussian'}

\item{GaussianDistrib}{list. Mean value and STD corresponding to the parameters sampled with gaussian distribution}

\item{ParmSet}{list. list of input parameters set to a specific value}

\item{nbSamples}{numeric. number of samples in training LUT}

\item{verbose}{boolean. when set to TRUE, prints message if hyperparameter adjustment performed during training}
}
\value{
InputPROSAIL
}
\description{
generate InputPROSAIL, following
\itemize{
\item either distribution defined in ATBD
\item or distribution defined in user
}
}
