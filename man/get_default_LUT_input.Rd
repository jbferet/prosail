% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_PROSAIL_LUT.R
\name{get_default_LUT_input}
\alias{get_default_LUT_input}
\title{This function sets default values for PROSAIL LUT simulation when not defined by user}
\usage{
get_default_LUT_input(
  TypeDistrib = NULL,
  GaussianDistrib = NULL,
  minval = NULL,
  maxval = NULL
)
}
\arguments{
\item{TypeDistrib}{list. specify if uniform or Gaussian distribution to be applied. default = Uniform}

\item{GaussianDistrib}{list. Mean value and STD corresponding to the parameters sampled with gaussian distribution}

\item{minval}{list. Defines the minimum value to be set for a list of parameters randomly produced}

\item{maxval}{list. Defines the maximum value to be set for a list of parameters randomly produced}
}
\value{
res list. list of default values corresponding to NULL input parameters
}
\description{
This function sets default values for PROSAIL LUT simulation when not defined by user
}
