---
title: '`prosail`: an R package to simulate canopy bi-directional reflectance 
factor from vegetation biophysical properties with the coupled model PROSAIL 
(PROSPECT + SAIL)'
tags:
  - R
  - vegetation radiative transfer modeling 
  - vegetation biophysical properties
  - remote sensing
  - prospect
  - sail
  - inversion
authors:
  - name: Jean-Baptiste Féret
    orcid: 0000-0002-0151-1334
    corresponding: true
    affiliation: 1
  - name: Florian de Boissieu
    orcid: 0000-0002-2185-9952
    affiliation: 1
affiliations:
 - name: TETIS, INRAE, AgroParisTech, CIRAD, CNRS, Université Montpellier, Montpellier, France
   index: 1
date: 26 March 2025
bibliography: paper.bib
---

# Summary

The PROSAIL model is a widely used radiative transfer model that combines the 
PROSPECT leaf optical properties model with the SAIL canopy bidirectional 
reflectance model [@jacquemoud2009; @berger2018].
It is designed to simulate the reflectance of vegetation canopies across a broad 
spectral domain, including the visible, near-infrared, shortwave infrared and 
thermal infrared regions of the electromagnetic spectrum.

PROSAIL integrates the optical properties of individual leaves, as calculated by 
PROSPECT [@jacquemoud1990], with the structural characteristics of the canopy, 
as modeled by SAIL [@verhoef1984]. 
This combination allows for the simulation of how light interacts with 
vegetation at both the leaf and canopy levels, taking into account factors such 
as leaf chemical content, leaf structure, canopy architecture, and viewing and 
illumination angles.

PROSAIL is commonly used in remote sensing applications to retrieve biophysical 
parameters of vegetation, such as leaf area index (LAI), leaf water and pigment 
content, the fraction of absorbed photosynthetically active radiation (fAPAR), 
albedo and fractional vegetation cover (fCover).
PROSAIL's ability to simulate bidirectional reflectance factor (BRF) under 
various conditions makes it a valuable tool for agricultural monitoring, 
ecological studies, and climate research.

Here, we introduce `prosail`, an R package which provides multiple versions of 
the model PROSAIL, coupled with the latest versions of PROSPECT distributed with 
the R package `prospect` [@feret2024]. 
The package includes functions to run simulations of BRF in forward mode for 
various optical sensors, as well as functions to apply user defined spectral 
response function (SRF), allowing simulation of any optical sensor with known or 
approximated BRF.
It also includes functions for model inversion using iterative 
optimization with and without prior information, as well as a module dedicated 
to hybrid inversion, combining physical modeling and machine learning 
regression, applicable to data tables and raster data. 

# Statement of need

The capacity to measure, map and monitor vegetation traits corresponding to 
biophysical and chemical properties is crucial to better understand ecosystem 
and agrosystem functions, as well as carbon, water and energy budgets. 
Vegetation radiative transfer models (RTMs) aim at describing the interactions 
between light and these biophysical and chemical properties, through absorption 
and scattering mechanisms. 
At leaf scale, the model PROSPECT (leaf optical PROperties SPECtra) [@feret2024] 
is currently the most popular physical model for the simulation of leaf optical 
properties and the accurate estimation of leaf chemistry 
[@feret2019; @spafford2021]. 
At canopy scale, various models can be employed, and the choice for a specific 
model may be driven by the complexity of the system to be simulated, as well as 
the capacity to describe the scene. 

The SAIL model is an example of four-stream representations of the radiative 
transfer equation, including two direct fluxes (incident solar flux and radiance 
in the viewing direction) and two diffuse fluxes (upward and downward 
hemispherical flux).
A system of four linear differential equations can then be analytically solved.
Detailed description of the functioning of the different SAIL versions can be 
found in [@verhoefbach2007] and [@verhoef2007].

More complex models include the Soil Canopy Observation of Photosynthesis and 
Energy fluxes (SCOPE) [@yang2021] model, which features layered description of 
the vegetation, enabling the simulation of vegetation with an understorey or 
with a vertical gradient in leaf properties, and the DART model [@gastellu2015], 
which allows for explicit 3D description of the canopy using geometric 
representations or point clouds derived from LiDAR acquisitions.

RTMs are increasingly used in remote sensing applications, in combination with 
machine learning: RTMs are used to produce a BRF dataset with corresponding 
vegetation properties, which are then used as training data to adjust a 
regression model with a machine learning algorithm. 
Such strategies, commonly referred to as hybrid methods [@verrelst2015], show 
multiple advantages : simulations are used in place of extensive in situ 
sampling, allowing regression model adjustment on virtually any source of 
optical imagery and any training sample size, while the use of machine learning 
ensures computation efficiency compared to iterative optimization when 
processing imagery data.

The relatively limited number of input parameters and the computational 
efficiency of the SAIL model makes it particularly interesting for operational 
applications, despite the main assumption of vegetation as a homogeneous turbid 
medium, where leaves are randomly distributed within the canopy, which is not 
accurate for row crops and heterogeneous canopies.

Various softwares currently allow application of hybrid inversion with PROSAIL 
simulations. 
The Sentinel Toolbox Application (SNAP) includes the `Biophysical Processor` 
module [@weiss2020], which provides an implementation of 
PROSAIL hybrid inversion combining PROSAIL simulations with an artificial neural 
network regression model, for the estimation of vegetation properties, such as 
LAI, fAPAR, fcover, canopy chlorophyll content and canopy water content. 
The sampling strategy for the definition of the training set, including the 
distribution and co-distribution of the input parameters, is described in the 
Algorithmic Theoretical Background Document (ATBD) [@weiss2020]) and cannot be 
adjusted by users.
Alternative distributions provide more interactive parameterization of the 
inversion strategy, such as the ARTMO toolbox [@rivera2014] developed in Matlab.

The R package `prosail` includes the versions of the model PROSPECT implemented 
in the package `prospect`: PROSPECT-D [@feret2017] and PROSPECT-PRO 
[@feret2021]. 
It also includes two versions of the model SAIL: 

- 4SAIL [@verhoef2007], a numerically stable adaptation of the SAILH version 
[@verhoef1998], which incorporate foliage hot spot to the original SAIL model

- 4SAIL2 [@verhoefbach2007], a two layers version of 4SAIL accommodating 
horizontal and vertical heterogeneities. 
This version includes additional parameters to account for crown clumping, and 
to describe the vertical distribution of two types of leaves described by 
specific leaf optical properties, which may describe different phenological 
stages or different vegetation types.

The R package `prosail` also provides a set of functions to help researchers 
aiming at experimenting for the estimation vegetation biophysical properties 
from optical remote sensing using RTM. 
These functions include : 

- the simulation of BRF for any optical sensor, based on their SRF

- the production of look-up-tables (LUTs) for simulated sensor BRF and 
corresponding PROSAIL input parameters, with or without additive and 
multiplicative noise

- the application of iterative optimization with possible added prior 
information on data tables

- the application of hybrid inversion on data tables and raster data


The package `prosail` does not intend to provide the same computational 
efficiency as SNAP. 
It does not provide a collection of models and methods as comprehensive as those 
provided with the ARTMO box neither. 
`prosail` provides a flexible framework to experiment with a hybrid inversion 
procedure with simple yet fast and efficient training stage, allowing 
experimenting on strategies for training data sampling, introduction of noise, 
feature selection over any type of optical sensor. 
It is particularly suitable for research and development, in order to identify 
potential and limitations of inversion strategies and parameterizations.
Resulting regression models can be applied on remote sensing data, but it may 
not be the appropriate software for scaling up vegetation monitoring 
applications. 

Alternative PROSAIL implementations can be found at 
[this webpage](http://teledetection.ipgp.jussieu.fr/prosail/).
This includes distributions in matlab, python and fortran programming languages. 
Note that alternative PROSAIL implementations are also available in packages 
written in [python](https://github.com/earth-chris/xleaf), 
[Julia](https://github.com/RemoteSensingTools/CanopyOptics.jl) and 
[R](https://github.com/ashiklom/rrtm).

# Overview

## PROSAIL simulation in forward mode

PROSAIL requires information intrinsic to vegetation : 

- leaf optical properties, including leaf directional hemispherical reflectance 
and transmittance. 
These leaf optical properties can be simulated with the R package `prospect`
[@feret2024] and readers are invited to refer to the documentation of this 
package for a comprehensive description of the leaf chemical and structure 
parameters accounted for by the various versions of the PROSPECT model. 

- Leaf Area Index (LAI), which is defined as the one-sided green leaf area per 
unit ground surface area in broadleaf canopies. 
It is dimensionless and is sometimes expressed in m<sup>2sup>.m<sup>-2sup>

- Leaf inclination distribution, which can be defined following different Leaf 
Inclination Distribution Functions (LIDF). 
The LIDF used in the original version of SAIL [verhoef1998] describes leaf
inclination based on two parameters: `LIDFa` controls the average leaf slope, 
and `LIDFb` controls the distribution's bimodality. 
Both parameters should be set between -1 and 1 and the sum of absolute values 
inferior to 1.
An alternative LIDF defined by [@campbell1990] can be applied, and requires a 
unique parameter corresponding to the average leaf angle.

- a foliage hot spot parameter which corresponds to the ratio of the correlation 
length of leaf projections in the horizontal plane and the canopy height 
[@kuusk1991; @breon2002; @verhoefbach2007]

PROSAIL also requires information extrinsic to vegetation : 

- the sun-observer geometry, defined by the sun zenith angle, the observer 
zenith angle, and the relative azimuth angle between sun and observer

- the soil reflectance

Additional parameters are required to run the version `4SAIL2` : 

- the optical properties corresponding to the second type of leaves 

- the fraction of LAI corresponding to a second type of leaves

- a layer dissociation factor between the main leaf layer and the second leaf 
layer

- the vertical crown cover percentage, which corresponds to the percentage of 
ground area covered with crowns as seen from nadir direction
 
- a tree shape factor corresponding to the ratio of crown diameter to crown 
height

SAIL produces a list of output variables corresponding to four top-of-canopy 
reflectance factors. 
Additional parameters required to compute albedo, fAPAR and fcover are also 
provided as outputs. 
The full list of outputs is then :

- `rddt`: bi-hemispherical reflectance factor

- `rsdt`: directional-hemispherical reflectance factor for solar incident flux

- `rdot`: hemispherical-directional reflectance factor in viewing direction

- `rsot`: bi-directional reflectance factor

- `abs_dir`: canopy absorptance for direct solar incident flux

- `abs_hem`: canopy absorptance for hemispherical diffuse incident flux

- `fcover`: Fraction of green Vegetation Cover (equals to 1 - beam transmittance 
in the target-view path)

- `rsdstar`: contribution of direct solar incident flux to albedo

- `rddstar`: contribution of hemispherical diffuse incident flux to albedo

The function `PRO4SAIL` using `4SAIL` can be called as follows : 

```r
# Load prosail package
library(prosail)

# define PROSPECT input variables. Refer to prospect tutorial for default values
input_prospect <- data.frame('CHL' = 40, 'CAR' = 8, 'ANT' = 0.0, 
                             'EWT' = 0.01, 'LMA' = 0.009, 'N' = 1.5)

# define input variables for SAIL. 
lai <- 5      # LAI
q <- 0.1      # foliage hot spot parameter
TypeLidf <- 2 # leaf inclination distribution function : Campbell
LIDFa <- 30   # average leaf angle (degrees)
tts <- 30     # geometry of acquisition: sun zenith angle (degrees)
tto <- 10     # geometry of acquisition: observer zenith angle (degrees)
psi <- 90     # geometry of acquisition: sun-observer azimuth (degrees)
rsoil <- SpecSOIL$Dry_Soil # soil reflectance

# run PROSAIL with 4SAIL
Ref_4SAIL <- PRO4SAIL(input_prospect = input_prospect, 
                      TypeLidf = TypeLidf, LIDFa = LIDFa, lai = lai,
                      q = q, tts = tts, tto = tto, psi = psi, rsoil = rsoil)
```

Additional input variables are required when using `4SAIL2` : 

```r
# define leaf chemical content corresponding to two types of leaves
input_prospect <- data.frame('CHL' = c(40,5), 'CAR' = c(8,5), 'ANT' = c(0,1), 
                             'EWT' = c(0.01,0.005), 'LMA' = c(0.009,0.008), 
                             'BROWN' = c(0.0,0.5), 'N' = c(1.5,2))

# define additional 4SAIL2 parameters
fraction_brown <- 0.5   # fraction of LAI corresponding to second type of leaves
diss <- 0.5             # layer dissociation factor between the leaf layers
Cv <- 1                 # vertical crown cover percentage
Zeta <- 1               # tree shape factor

# run PROSAIL with 4SAIL2
Ref_4SAIL2 <- PRO4SAIL(SAILversion = '4SAIL2', input_prospect = input_prospect,
                       TypeLidf = TypeLidf, LIDFa = LIDFa, lai = lai,
                       q = q, tts = tts, tto = tto, psi = psi, rsoil = rsoil,
                       fraction_brown = fraction_brown, diss = diss, 
                       Cv = Cv, Zeta = Zeta)
```

## Computation of BRF from PROSAIL outputs

The BRF results from the combination of hemispherical-directional and 
bi-directional reflectance factors with the relative contribution of direct and
diffuse radiation. 
The direct and diffuse radiations are defined as described by [@francois2002]. 
Users can then compute BRF by setting a global contribution of diffuse skylight 
radiation with the `skyl` parameter, which is the share of diffuse flux in 
global radiation. 
Alternatively, the direct and diffuse radiations can taken into account by 
computing `skyl` based on the equation proposed by [@spitters1986], knowing the 
sun zenith angle and assuming clear sky conditions. 

```r
# Compute BRF with known skyl
BRF_4SAIL <- Compute_BRF(rdot = Ref_4SAIL$rdot, rsot = Ref_4SAIL$rsot,
                         skyl = 0.23, SpecATM_Sensor = SpecATM)

# Compute BRF assuming clear sky conditions and using sun zenith angle
BRF_4SAIL <- Compute_BRF(rdot = Ref_4SAIL$rdot, rsot = Ref_4SAIL$rsot,
                         tts = 40, SpecATM_Sensor = SpecATM)
```

## Computation of fAPAR and albedo

The fAPAR in the direction of the observer can be computed by combining both 
canopy absorptance for direct solar incident flux and canopy absorptance for 
hemispherical diffuse incident flux.

```r
fAPAR_4SAIL <- Compute_fAPAR(abs_dir = Ref_4SAIL$abs_dir,
                             abs_hem = Ref_4SAIL$abs_hem,
                             tts = tts, SpecATM_Sensor = SpecATM)
```

Finally, the albedo can be derived from direct solar incident flux and 
hemispherical diffuse incident flux.

```r
albedo_4SAIL <- Compute_albedo(rsdstar = Ref_4SAIL$rsdstar,
                               rddstar = Ref_4SAIL$rddstar,
                               tts = tts, SpecATM_Sensor = SpecATM)
```

## Simulating sensor BRF

The simulation of sensor BRF requires the SRF of the sensor. 
A selection of SRF corresponding to multiple sensors is already implemented in 
`prosail` :

- Sentinel-2A, 2B and 2C

- Landsat-7, 8 and 9

- MODIS

- Venus

- SPOT-6 / 7

- Pleiades 1

Users can also define the sensor of their choice, either by providing central 
wavelength and full width half maximum (fwhm) corresponding to each band and 
assuming Gaussian response for each band, or by providing a file describing the 
exact SRF. 

```r
# get the spectral response function (SRF) for Sentinel-2A
SRF_S2 <- GetRadiometry('Sentinel_2A')

# get the SRF corresponding to a sensor defined by user. 
# Hyperspectral sensor: 10 nm spectral sampling & 10 nm fwhm for all bands
wl <- seq(400, 2500, by = 10)
fwhm <- rep(x = 10, length = length(wl))
spectral_properties = data.frame(wl = wl, fwhm = fwhm)
sensor_name <- 'Hyperspectral_Sensor'
SRF_HSI <- prosail::GetRadiometry(sensor_name = sensor_name,
                                  spectral_properties = spectral_properties)
```

## PROSAIL inversion through iterative optimization

### Iterative optimization with no prior information

PROSAIL is a relatively simple and computationally efficient model. 
As with PROSPECT, inversion based on iterative optimization can be considered. 
However, unlike PROSPECT inversion which is well-posed in most situations, 
PROSAIL inversion is inherently ill-posed. 
Therefore, standard iterative optimization algorithms consisting in minimizing a 
cost function will converge towards local minimum, unless global optimization 
algorithms are applied. 
This means that initial conditions of the iterative optimization influence the 
output of a local optimization.

Multiple inversion strategies based on iterative optimization are implemented in 
the package `prospect`. 
The main iterative optimization algorithm is based on the minimization of a 
multivariable function with nonlinear constraints. 
This procedure is based on the function `fmincon` included in the package 
`pracma`.
The default criterion for the cost function corresponds to the root mean square 
error (RMSE) between a reference BRF and BRF simulated with PROSAIL. 
The example code below provides an example of iterative optimization applied to
simulated BRF at native spectral sampling (1 nm), and to Sentinel-2 BRF. 

```r
# define initial value, lower and upper bounds for inversion
Init <- data.frame('CHL' = 40, 'CAR' = 10, 'EWT' = 0.01, 'LMA' = 0.01, 
                   'lai' = 3, 'LIDFa' = 50, 'N' = 1.5)
LB <- data.frame('CHL'=5, 'CAR'=1, 'EWT'=0.002, 'LMA'= 0, 
                 'lai' = 0.5, 'LIDFa' = 30, 'N' = 1)
UB <- data.frame('CHL' = 80, 'CAR' = 20, 'EWT' = 0.03, 'LMA' = 0.03, 
                 'lai' = 6, 'LIDFa' = 80, 'N' = 3)

# define parameters set for inversion
parm_set <- data.frame('tts' = 40, 'tto' = 0, 'psi' = 60,  'psoil' = 0,
                      'ANT' = 0, 'BROWN' = 0, 'q' = 0.1)

# compute soil reflectance
rsoil <- parm_set$psoil*SpecSOIL$Dry_Soil+(1-parm_set$psoil)*SpecSOIL$Wet_Soil

# simulate canopy BRF with 1 nm sampling
truth <- data.frame('CHL' = 60, 'CAR' = 8, 'EWT' = 0.015, 'LMA' = 0.005,
                    'lai' = 3,  'LIDFa' = 60, 'N' = 1.8)
Refl_1nm <- PRO4SAIL(N = truth$N, CHL = truth$CHL, CAR = truth$CAR,
                     ANT = parm_set$ANT, BROWN = parm_set$BROWN, EWT = truth$EWT,
                     LMA = truth$LMA, TypeLidf = 2, lai = truth$lai,
                     q = parm_set$q, LIDFa = Init$LIDFa, rsoil = rsoil,
                     tts = parm_set$tts, tto = parm_set$tto, psi = parm_set$psi)
brf_1nm <- prosail::Compute_BRF(rdot = Refl_1nm$rdot, rsot = Refl_1nm$rsot,
                                tts = parm_set$tts, SpecATM_Sensor = SpecATM)

# invert PROSAIL on BRF with 1 nm spectral sampling
est_1nm <- Invert_PROSAIL(brf_mes = brf_1nm$BRF, initialization = Init,
                          lower_bound = LB, upper_bound = UB,
                          SpecPROSPECT_Sensor = SpecPROSPECT_FullRange,
                          SpecATM_Sensor = SpecATM, SpecSOIL_Sensor = SpecSOIL,
                          TypeLidf = 2, parm_set = parm_set)

# convert spectral input based on Sentinel-2 SRF
SRF_S2 <- GetRadiometry('Sentinel_2A')
brf_S2 <- applySensorCharacteristics(wvl = prospect::SpecPROSPECT_FullRange$lambda, 
                                     SRF = SRF_S2, InRefl = brf_1nm$BRF)
SpecPROSPECT_s2 <- applySensorCharacteristics(wvl = prospect::SpecPROSPECT_FullRange$lambda, 
                                           SRF = SRF_S2, InRefl = prospect::SpecPROSPECT_FullRange)
SpecSOIL_s2 <- applySensorCharacteristics(wvl = prospect::SpecPROSPECT_FullRange$lambda, 
                                           SRF = SRF_S2, InRefl = SpecSOIL)
SpecATM_s2 <- applySensorCharacteristics(wvl = prospect::SpecPROSPECT_FullRange$lambda, 
                                           SRF = SRF_S2, InRefl = SpecATM)

# invert PROSAIL on Sentinel-2 BRF
est_s2 <- Invert_PROSAIL(brf_mes = brf_S2, initialization = Init,
                         lower_bound = LB, upper_bound = UB,
                         SpecPROSPECT_Sensor = SpecPROSPECT_s2,
                         SpecATM_Sensor = SpecATM_s2, SpecSOIL_Sensor = SpecSOIL_s2,
                         TypeLidf = 2, parm_set = parm_set)
```

### Iterative optimization with prior information

Prior information can be added to the cost function in order to regularize the 
ill-posed inverse problem [@combal2003].
Multiple possibilities exist to account for prior information: 

- the definition of `initialization` variable according to initial assumptions

- the addition of a term corresponding to the probability density function of 
part or all of the parameters to be estimated in the cost function. 

The code blow illustrates the introduction of prior information on LIDF, with the 
definition of a probability density function corresponding to the average leaf 
angle `LIDFa`. 
A weighting factor associated to the relative importance of this prior 
information with respect to the main criterion to minimize (RMSE between 
reference BRF and simulated BRF) can also be adjusted.


```r
# define prior information on average leaf angle LIDFa
prior_info <- list('mean' = data.frame('LIDFa' = 60), 
				   'sd' = data.frame('LIDFa' = 10), 
				   'WeightPrior' = 0.01)


# invert PROSAIL on BRF with 1 nm spectral sampling and probability density 
# function for LIDFa
est_1nm_p <- Invert_PROSAIL(brf_mes = brf_1nm$BRF, initialization = Init,
                            lower_bound = LB, upper_bound = UB,
                            SpecPROSPECT_Sensor = SpecPROSPECT_FullRange,
                            SpecATM_Sensor = SpecATM, 
                            SpecSOIL_Sensor = SpecSOIL,
                            TypeLidf = 2, parm_set = parm_set,
                            prior_info =  prior_info)

# invert PROSAIL on Sentinel-2 BRF and probability density function for LIDFa
est_s2_p <- Invert_PROSAIL(brf_mes = brf_S2, initialization = Init,
                           lower_bound = LB, upper_bound = UB,
                           SpecPROSPECT_Sensor = SpecPROSPECT_s2,
                           SpecATM_Sensor = SpecATM_s2, 
                           SpecSOIL_Sensor = SpecSOIL_s2,
                           TypeLidf = 2, parm_set = parm_set,
                           prior_info =  prior_info)
```

## PROSAIL hybrid inversion

Hybrid model inversion is defined here as a two-step inversion procedure: 

- run RTM to produce a simulated LUT including sensor BRF and corresponding 
PROSAIL input parameters, including vegetation properties. 

- train a machine learning regression algorithm to estimate a vegetation 
property or a set of vegetation properties of interest from sensor BRF.

The next subsections describe the strategy followed by the hybrid inversion 
implemented in `prosail`.

### Simulating BRF data to prepare for regression model training 

Several factors need to be accounted for when preparing for the training stage 
of the regression models : 

- the sampling strategy (range, distribution, co-distributions for input 
parameters)

- the addition of additive and multiplicative noise to the simulated BRF

- the definition of spectral bands to be used for the training of the regression
model

The default parameterization of `prosail` hybrid inversion relies on the 
production of a BRF LUT compliant with the specifications defined in the ATBD 
for the Biophysical Processor of the Sentinel toolbox [@weiss2020]. 
However, each step of the simulation of a training LUT can be defined by user.


### Training a machine learning regression algorithm

The machine learning strategy implemented in `prosail` consists in a 
parsimonious ensemble method based on a bootstrap aggregating (bagging) 
prediction of biophysical properties from multiple support vector regression 
(SVR) models.
Each SVR model is trained with a limited number of samples, ensuring fast 
training stage. 
The predicted value then corresponds to the mean prediction from a set of 
SVR models.
The standard deviation can also be derived from this ensemble of predictors.
However, standard deviation may not provide accurate model uncertainty 
quantification [@palmer2022].

The default parameterization of this inversion strategy uses an ensemble of 20 
SVR models, each of them trained with 100 samples. 
Therefore, only 2000 samples are required for the training stage, which is 
relatively low compared to the 41472 cases used for the training of the 
artificial neural networks used in the SNAP toolbox. 
This parsimonious ensemble method ensures very fast training stage and similar 
performances as those obtained with the SNAP toolbox. 

### Performing full hybrid inversion

Figure \ref{fig:hybrid} summarizes the workflow applied to perform hybrid 
inversion. 
The function `train_prosail_inversion` combines all aforementioned steps to 
produce BRF simulations and  train SVR models. 
It provides regression models as outputs, which can then be directly used to 
estimate vegetation biophysical properties from data tables or raster data. 

![Workflow of PROSAIL hybrid inversion implemented in the package prosail. \label{fig:hybrid}](FlowChart_HybridInversion_JOSS.tif){ width=75% }

The example code below illustrates the training stage with the function
`train_prosail_inversion`, to estimate LAI, fCover and fAPAR from Sentinel-2, 
using the three spectral bands corresponding to green channel (B3), red channel 
(B4) and near infrared channel (B8) as specified in the ATBD document.

```r
# define sensor to simulate with PROSAIL
SRF <- get_radiometry('Sentinel_2')

# define parameters to estimate
Parms2Estimate <- c('lai', 'fCover', 'fAPAR')

# define spectral bands required to train SVR model for each variable
Bands2Select <- list()
S2BandSelect <- c('B3','B4','B8')
for (parm in Parms2Estimate)
  Bands2Select[[parm]] <- match(S2BandSelect,SRF$Spectral_Bands)

# define output directory where LUTs will be saved
PROSAIL_ResPath <- './HybridInversion'
dir.create(path = PROSAIL_ResPath, showWarnings = FALSE,recursive = TRUE)

# define ranges for geometry of acquisition
GeomAcq <- list('min' = data.frame('tto' = 0, 'tts' = 20, 'psi' = 0), 
                'max' = data.frame('tto' = 10, 'tts' = 30, 'psi' = 360))

# train model
modelSVR <- train_prosail_inversion(Parms2Estimate = Parms2Estimate,
                                    atbd = TRUE, GeomAcq = GeomAcq, 
                                    SRF = SRF, 
                                    Bands2Select = Bands2Select, 
                                    output_dir = PROSAIL_ResPath)
```

The processing steps called in `train_prosail_inversion` can be broken down as 
follows:


```r
# produce input variables following ATBD specifications and geometry of acq.
input_prosail <- get_input_prosail(atbd = TRUE, GeomAcq = GeomAcq)

# generate a LUT with 1 nm spectral sampling
res <- Generate_LUT_PROSAIL(SAILversion = '4SAIL', input_prosail = input_prosail,
                            SpecPROSPECT = prospect::SpecPROSPECT_FullRange,
                            SpecSOIL = SpecSOIL, SpecATM = SpecATM)

# include fcover and fAPAR (derived from SAIL reflectance / absortion factors)
# in the pool of variables to explain
input_prosail$fCover <- res$fCover
input_prosail$fAPAR <- res$fAPAR

# apply Sentinel-2 SRF to get sensor BRF
BRF_LUT <- applySensorCharacteristics(wvl = prospect::SpecPROSPECT_FullRange$lambda, 
                                      SRF = SRF, InRefl = res$BRF)
# identify spectral bands in LUT
rownames(BRF_LUT) <- SRF$Spectral_Bands

# apply noise and produce a LUT for each parameter (may vary with variable to explain)
for (parm in Parms2Estimate) 
  BRF_LUT_Noise[[parm]] <- apply_noise_atbd(BRF_LUT)

# train a set of SVR models for each parameter
modelSVR <- list()
for (parm in Parms2Estimate)
  modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                           input_variables = input_prosail[[parm]])
```

The resulting models can then be applied on data tables or raster data. 
The following example perform predictions based on models trained previously. 
The prediction returns mean value obtained form the ensemble of regression 
models for each sample, as well as corresponding standard deviation.

```r
MeanEstimate <- StdEstimate <- list()
for (parm in Parms2Estimate){
  HybridRes <- PROSAIL_Hybrid_Apply(RegressionModels = modelSVR[[parm]],
                                    Refl = BRF_LUT_Noise[[parm]])
  MeanEstimate[[parm]] <- HybridRes$MeanEstimate
  StdEstimate[[parm]] <- HybridRes$StdEstimate
}
```


### Application to imagery data

Here, we illustrate how to run `prosail` hybrid inversion on Sentinel-2 imagery.
This requires additional packages, including the package `sf` to handle vector 
data, and the package [`preprocs2`](https://jbferet.gitlab.io/preprocs2/), which 
is dedicated to downloading and preprocessing Sentinel-2 data. 

Sentinel-2 data can be accessed using various modalities and providers.
`preprocs2` allows access to Sentinel-2 data from STAC API. 
However, for the sake of comparison between biophysical variables produced from 
the SNAP toolbox and those produced with the `prosail` hybrid inversion, the 
Sentinel-2 Level-2A product corresponding to tile 30SWJ acquired on May 13th, 
2021 was downloaded from the 
[Copernicus browser](https://browser.dataspace.copernicus.eu).

Once the SAFE product downloaded, `preprocs2` is used to crop and save 
reflectance data corresponding to an area of interest defined within the 
Sentinel-2 tile footprint. 


```r
# libraries
library(prosail)
library(preprocS2)
library(sf)

# define path for original image and outputs
main_dir <- './barrax'                              # main directory
safe_path <- file.path(main_dir,                    # path to SAFE L2A S2 image
                       'S2B_MSIL2A_20210513T105619_N0500_R094_T30SWJ_20230228T104126.SAFE')
out_s2 <- file.path(main_dir, 'S2_subset')          # S2 subset directory
out_BP <- file.path(main_dir, 'prosail_inversion')  # prosail products directory
dir.create(path = out_BP, showWarnings = FALSE,recursive = TRUE)

# define area of interest included in S2 tile
aoi_bbox <- st_bbox(obj = c('xmin' = 571626, 'ymin' = 4324941, 
                            'xmax' = 582995, 'ymax' = 4331253))
aoi <- bbox_to_poly(aoi_bbox, crs = 32630)
path_aoi <- file.path(main_dir, 'barrax.gpkg')
st_write(obj = aoi, dsn = path_aoi, driver = 'GPKG', delete_dsn = T)

s2_products <- extract_from_safe(safe_path = safe_path, 
                                 path_aoi = path_aoi, 
                                 output_dir = out_s2)
```


### Comparison with SNAP toolbox

LAI, fcover and fAPAR computed from Sentinel-2 images with the hybrid inversion 
implemented in `prosail` are compared to those produced with SNAP in Figure 
\ref{fig:SNAP_prosail}. 

![Biophysical properties estimated with `prosail` vs Biophysical properties estimated with SNAP. \label{fig:SNAP_prosail}](compare_SNAP_prosail.png){ width=90% }

The two methods show relatively good consistency with Pearson correlation 
coefficient ranging between 0.98 and 0.99 for the three variables.
Differences remain between the two implementations, including differences in the 
version of the PROSPECT model (SNAP uses a version anterior to PROSPECT-D), and 
differences in the soil reflectance used to produce the BRF LUT. 
These differences between the two approaches contribute to explain the 
differences observed between corresponding biophysical properties, particularly 
for lower vegetation cover, for which BRF is more influenced by soil reflectance. 

# Conclusion

We introduce `prosail`, an R package dedicated to the canopy reflectance model 
PROSAIL, coupling PROSPECT and SAIL. 
`prosail` is coupled with the R package `prospect` in order to allow 
seamless integration of future versions of the leaf model. 
It allows simulation of bi-directional reflectance factor for any type of 
optical sensor if the spectral response function is available, including 
multispectral sensors and hyperspectral sensors. 
The package also includes a collection of inversion procedures, including 
iterative optimization with and without prior information, and hybrid inversion 
applicable to data tables and raster data. 

The estimation of vegetation biophysical properties with `prosail` hybrid 
inversion is consistent with estimations produced with the method implemented in
the Sentinel toolbox, which is currently the reference method. 

`prosail` hybrid inversion does not intend to be computationally as efficient as 
SNAP and may not meet requirements for large scale vegetation monitoring 
applications. 
`prosail` hybrid inversion offers a fully adjustable hybrid inversion framework, 
including an original parsimonious machine learning regression method, allowing 
fast and efficient training for any type of sensor. 
It is a valuable tool for experimenting on the potential and limitation of 
physically-based retrieval of vegetation traits from various types of optical 
sensors, including airborne imaging spectroscopy, as well satellite missions 
already operational and in preparation. 

# Availability

`prosail` is an open-source software package made available under the MIT license. 
Tutorials are available at [https://jbferet.gitlab.io/prosail/](https://jbferet.gitlab.io/prosail/).

# Acknowledgements

The authors acknowledge financial support from Agence Nationale de la Recherche 
(BioCop project — ANR-17-CE32-0001).
We are grateful to Wout Verhoef for the development of the initial version of 
the 4SAIL and 4SAIL2 models. 

# References
