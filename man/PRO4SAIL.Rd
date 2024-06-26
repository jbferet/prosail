% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_PROSAIL.R
\name{PRO4SAIL}
\alias{PRO4SAIL}
\title{PROSAIL simulation based on a set of combinations of input parameters}
\usage{
PRO4SAIL(
  Spec_Sensor = NULL,
  Input_PROSPECT = NULL,
  N = 1.5,
  CHL = 40,
  CAR = 8,
  ANT = 0,
  BROWN = 0,
  EWT = 0.01,
  LMA = NULL,
  PROT = 0,
  CBC = 0,
  alpha = 40,
  TypeLidf = 2,
  LIDFa = 60,
  LIDFb = NULL,
  lai = 3,
  q = 0.1,
  tts = 30,
  tto = 0,
  psi = 60,
  rsoil = NULL,
  fraction_brown = 0,
  diss = 0,
  Cv = 1,
  Zeta = 1,
  SAILversion = "4SAIL",
  BrownLOP = NULL
)
}
\arguments{
\item{Spec_Sensor}{list. Includes optical constants required for PROSPECT
refractive index, specific absorption coefficients and spectral bands}

\item{Input_PROSPECT}{list. PROSPECT input variables}

\item{N}{numeric. Leaf structure parameter}

\item{CHL}{numeric. Chlorophyll content (microg.cm-2)}

\item{CAR}{numeric. Carotenoid content (microg.cm-2)}

\item{ANT}{numeric. Anthocyain content (microg.cm-2)}

\item{BROWN}{numeric. Brown pigment content (Arbitrary units)}

\item{EWT}{numeric. Equivalent Water Thickness (g.cm-2)}

\item{LMA}{numeric. Leaf Mass per Area (g.cm-2)}

\item{PROT}{numeric. protein content  (g.cm-2)}

\item{CBC}{numeric. NonProtCarbon-based constituent content (g.cm-2)}

\item{alpha}{numeric. Solid angle for incident light at surface of leaf}

\item{TypeLidf}{numeric. Type of leaf inclination distribution function}

\item{LIDFa}{numeric.
if TypeLidf ==1, controls the average leaf slope
if TypeLidf ==2, controls the average leaf angle}

\item{LIDFb}{numeric.
if TypeLidf ==1, controls the distribution's bimodality
if TypeLidf ==2, unused}

\item{lai}{numeric. Leaf Area Index}

\item{q}{numeric. Hot Spot parameter}

\item{tts}{numeric. Sun zeith angle}

\item{tto}{numeric. Observer zeith angle}

\item{psi}{numeric. Azimuth Sun / Observer}

\item{rsoil}{numeric. Soil reflectance}

\item{fraction_brown}{numeric. Fraction of brown leaf area}

\item{diss}{numeric. Layer dissociation factor}

\item{Cv}{numeric. vertical crown cover percentage
= \% ground area covered with crowns as seen from nadir direction}

\item{Zeta}{numeric. Tree shape factor
= ratio of crown diameter to crown height}

\item{SAILversion}{character. choose between 4SAIL and 4SAIL2}

\item{BrownLOP}{dataframe. optical properties for brown vegetation
\itemize{
\item WVL, Reflectance, Transmittance
\item Set to NULL if use PROSPECT to generate it
}}
}
\value{
list. rdot,rsot,rddt,rsdt
rdot: hemispherical-directional reflectance factor in viewing direction
rsot: bi-directional reflectance factor
rsdt: directional-hemispherical reflectance factor for solar incident flux
rddt: bi-hemispherical reflectance factor
}
\description{
PROSAIL simulation based on a set of combinations of input parameters
}
