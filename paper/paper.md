---
title: '`prosail`: an R package to simulate canopy reflectance with the coupled 
model PROSAIL (PROSPECT + SAIL)'
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
date: 25 July 2025
bibliography: paper.bib
---

# Summary

PROSAIL is a widely used vegetation radiative transfer model (RTM), coupling the 
*leaf optical PROperties SPECtra* (PROSPECT) model, which simulates leaf optical 
properties [@jacquemoud1990], with the *Scattering by Arbitrarily Inclined Leaves* 
model (SAIL) which simulates canopy bidirectional reflectance 
[@verhoef1984; @jacquemoud2009; @berger2018] across the visible to shortwave 
infrared regions of the electromagnetic spectrum. 
PROSAIL's ability to simulate plant canopy spectral and directional reflectance 
under various conditions makes it a valuable tool for agriculture, ecology, and 
climate research [@poulter2023]. 
It is used in combination with inversion procedures to retrieve vegetation 
biophysical properties, including leaf area index (LAI), leaf water and 
pigment content, the fraction of absorbed photosynthetically active radiation 
(fAPAR), and fractional vegetation cover (fCover) among others.
The capacity to retrieve a biophysical property mainly depends on the spectral 
characteristics of the sensor, and the acquisition of sufficient spectral 
information to solve an ill-posed problem. 

We introduce `prosail`, an R package which provides multiple versions of 
the model PROSAIL, coupled with the latest versions of PROSPECT distributed with 
the R package `prospect` [@feret2024]. 
`prosail` can simulate surface reflectance in forward mode for various optical 
sensors, based on their spectral response function (SRF).
It also includes functions for model inversion using iterative 
optimization, and a module dedicated to hybrid inversion, combining physical 
modeling and ML regression. 

# Statement of need

The capacity to measure, map and monitor vegetation traits is crucial to better 
understand ecosystem and agrosystem functions, carbon, water and energy budgets. 
Vegetation RTMs aim at describing the interactions between light and vegetation
based on their biophysical and chemical properties, through absorption and 
scattering mechanisms. 
At leaf scale, PROSPECT [@feret2024] is currently the most popular physical 
model for the simulation of leaf optical properties and the accurate estimation 
of leaf chemistry [@feret2019; @spafford2021]. 
At canopy scale, various models can be employed, and the choice for a specific 
model may be driven by the complexity of the system to be simulated, and
the capacity to describe the scene. 

The SAIL model is an example of four-stream representations of the radiative 
transfer equation, including two direct fluxes (incident solar flux and radiance 
in the viewing direction) and two diffuse fluxes (upward and downward 
hemispherical flux).
A system of four linear differential equations can then be analytically solved.
Detailed description of the functioning of the different SAIL versions can be 
found in [@verhoefbach2007] and [@verhoef2007].

More complex models include the *Soil Canopy Observation of Photosynthesis and 
Energy fluxes* (SCOPE) [@yang2021] model, which features layered description of 
the vegetation, enabling the simulation of vegetation with an understorey or 
with a vertical gradient in leaf properties, and the 
*Discrete Anisotropic Radiative Transfer* DART model [@gastellu2015], which 
allows for explicit 3-dimensional description of the canopy using geometric 
representations with triangular meshes or turbid voxels that can be derived from 
LiDAR acquisitions.

RTMs are increasingly used in remote sensing applications, in combination with 
machine learning (ML): RTMs produce a surface reflectance dataset with 
corresponding vegetation properties, which are then used as training data to 
adjust a regression model with a ML algorithm. 
These hybrid strategies show multiple advantages [@verrelst2015]: simulations 
are used in place of extensive in situ sampling, while the use of ML 
ensures computation efficiency compared to iterative optimization.

The relatively limited number of input parameters and the computational 
efficiency of the SAIL model is interesting for operational applications, 
despite the main assumption of vegetation as a homogeneous turbid medium, where 
leaves are randomly distributed within the canopy, which is not accurate for 
row crops and heterogeneous canopies.

# State of the field

Various softwares currently provide procedures for hybrid inversion with PROSAIL 
simulations. 
The *Sentinel Toolbox Application Platform* (SNAP) includes the 
`Biophysical Processor` module [@weiss2020], with PROSAIL hybrid inversion 
combining PROSAIL simulations with an artificial neural network regression 
model, for the estimation of vegetation properties, including LAI, fAPAR, 
fcover, canopy chlorophyll and water content. 
The sampling of the training set, including the distribution and co-distribution 
of the input parameters, is described in the 
*Algorithmic Theoretical Background Document* (ATBD) [@weiss2020]).
Alternative distributions provide more interactive parameterization of the 
inversion strategy, such as the *Automated Radiative Transfer Models Operator*
(ARTMO) toolbox [@rivera2014] developed in Matlab.

The R package `prosail` includes the versions of the model PROSPECT implemented 
in the package `prospect`: PROSPECT-D [@feret2017] and PROSPECT-PRO 
[@feret2021]. 
It also includes two versions of the model SAIL, *4SAIL* [@verhoef2007], and 
*4SAIL2* [@verhoefbach2007], a two layers version of *4SAIL*. 

The package `prosail` does not intend to provide the same computational 
efficiency as SNAP. 
It does not provide a collection of models and methods as comprehensive as those 
provided with the ARTMO box neither. 
`prosail` provides a flexible open source framework to experiment with hybrid 
inversion procedures, using simple yet fast and efficient training stage, 
allowing experimenting on strategies for training data sampling, introduction 
of noise, feature selection over any type of optical sensor. 
It is particularly suitable for research and development, in order to identify 
potential and limitations of inversion strategies and parameterizations.
Resulting regression models can be applied on remote sensing data, but it may 
not be the appropriate software for scaling up vegetation monitoring 
applications. 

Alternative PROSAIL implementations can be found at 
[this webpage](http://teledetection.ipgp.jussieu.fr/prosail/).
This includes distributions in matlab, python and fortran programming languages. 
Note that alternative PROSAIL implementations are also available in packages 
written in [python](https://github.com/jgomezdans/prosail), 
[Julia](https://github.com/RemoteSensingTools/CanopyOptics.jl) and 
[R](https://github.com/ashiklom/rrtm).

# Software design

`prosail` uses the R package `prospect` [@feret2024] for the 
simulation of leaf optical properties. 
It provides a user-friendly and modular set of functions, to run 
simulations of vegetation canopy reflectance, and to predict vegetation 
biophysical properties using various inversion strategies. 
Model inversion is a multi-step procedure requiring simulation of sensor 
reflectance to train a ML regression algorithm, then 
applicable to any new data source, including raster data from airborne and 
spaceborne sensors. 

- The generation of realistic reflectance simulations requires definition of 
input parameter sampling strategy, accounting for sensor SRF and uncertainties. 
`prosail` provides functions to produce this training dataset as a unique step 
following standard parameterization, or step by step, in order for users to 
understand the logic, or adjust each step to their needs.

- A limited number of ML algorithms is currently provided, but the 
software design will allow integration of additional algorithms in future 
versions.

- Raster processing is mainly handled with the `terra` [@terra] package, to 
ensure compatibility with most raster data formats. 

- `prosail` can be used in combination with 
[`preprocS2`](https://jbferet.gitlab.io/preprocs2/), a package dedicated to 
downloading and preprocessing of Earth observation data following the 
[STAC](https://stacspec.org/en) specification, in order to produce a fully 
automated image access and processing workflow.

# Research impact statement

`prosail` has been used in multiple research publications since its early 
developments [@hauser2021], [@ferreira2026], [@feret2026].
It is currently used in multiple research projects, and is actively maintained. 


# Overview

## PROSAIL simulation in forward mode

PROSAIL requires information intrinsic to vegetation : 

- leaf directional hemispherical reflectance and transmittance. 
Readers are invited to refer to the documentation of `prospect` for 
a comprehensive description of the leaf chemical and structure 
parameters accounted for. 
Measured leaf optical properties can also be used instead of simulated ones.

- LAI, which is the one-sided green leaf area per unit ground surface area in 
broadleaf canopies. 
It is dimensionless, sometimes expressed in $m^{2}.m^{-2}$

- Leaf inclination distribution, which can be defined following different 
*Leaf Inclination Distribution Functions* (LIDF). 

- a foliage hot spot parameter which corresponds to the ratio of the correlation 
length of leaf projections in the horizontal plane and the canopy height 
[@kuusk1991; @breon2002; @verhoefbach2007]

PROSAIL also requires information extrinsic to vegetation : 

- the sun-observer geometry, defined by the sun zenith angle, the observer 
zenith angle, and the relative azimuth angle between sun and observer

- the soil reflectance

SAIL produces four top-of-canopy reflectance factors and additional parameters, 
required to compute albedo, fAPAR and fcover :

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

The function `prosail` using `4SAIL` can be called as follows : 

```r
# Load prosail package
library(prosail)

# define PROSPECT input variables
input_prospect <- data.frame('chl' = 40, 'car' = 8, 'ant' = 0.0, 
                             'ewt' = 0.01, 'lma' = 0.009, 'n_struct' = 1.5)

# define input variables for SAIL. 
lai <- 5        # LAI
hotspot <- 0.1  # foliage hot spot parameter
type_lidf <- 2  # leaf inclination distribution function : Campbell
lidf_a <- 30    # average leaf angle (degrees)
tts <- 30       # sun zenith angle (degrees)
tto <- 10       # observer zenith angle (degrees)
psi <- 90       # sun-observer azimuth (degrees)
rsoil <- spec_soil_ossl $soil_01 # soil reflectance selected from OSSL library

# run PROSAIL with 4SAIL
ref_4sail <- prosail(input_prospect = input_prospect, 
                      type_lidf = type_lidf, lidf_a = lidf_a, lai = lai,
                      hotspot = hotspot, tts = tts, tto = tto, psi = psi, 
                      rsoil = rsoil)
                      
```

## Simulating surface reflectance acquired by a sensor 

The simulation of sensor surface reflectance requires the SRF of the sensor. 
A selection of SRF corresponding to multiple sensors is already implemented in 
`prosail` :

- Sentinel-2A, 2B and 2C

- Landsat-7, 8 and 9

- MODIS

- Venus

- SPOT-6 / 7

- Pleiades 1

Users can also define the sensor of their choice, either by providing central 
wavelength and full width at half maximum (fwhm) corresponding to each band and 
assuming gaussian response for each band, or by providing a file describing the 
exact SRF. 


## PROSAIL hybrid inversion

Hybrid model inversion is a two-step inversion procedure: 

- run RTM to produce a simulated LUT including sensor surface reflectance and 
corresponding PROSAIL input parameters, including vegetation properties. 

- train a ML regression algorithm to estimate a vegetation 
property or a set of vegetation properties of interest from sensor reflectance.

The next subsections describe the strategy followed by the hybrid inversion 
implemented in `prosail`.

### Simulating surface reflectance to prepare for regression model training 

Several factors need to be accounted for when preparing for the training stage 
of the regression models : 

- the sampling strategy (range, distribution, co-distributions for input 
parameters)

- the addition of noise to the simulated surface reflectance

- the definition of spectral bands to be used for the training of the regression
model

The default parameterization of `prosail` hybrid inversion relies on the 
production of a surface reflectance LUT compliant with the specifications 
defined in the ATBD for the Biophysical Processor of the Sentinel toolbox [@weiss2020]. 
However, each step of the simulation of a training LUT can be defined by user.


### Training a ML regression algorithm

The ML strategy implemented in `prosail` consists in a 
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

The default parameterization of this inversion strategy uses an ensemble of 10 
SVR models, each of them trained with 200 samples. 
Therefore, only 2000 samples are required for the training stage, which is 
relatively low compared to the 41472 cases used for the training of the 
artificial neural networks used in the SNAP toolbox. 
This parsimonious ensemble method ensures very fast training stage and similar 
retrieval performances as those obtained with the SNAP toolbox. 

### Performing full hybrid inversion

Figure \ref{fig:hybrid} summarizes the workflow applied to perform hybrid 
inversion. 
The function `train_prosail_inversion` combines all aforementioned steps to 
produce surface reflectance simulations and train SVR models. 
It produces regression models which can then be directly used to estimate 
vegetation biophysical properties from data tables or raster data. 

![Workflow of PROSAIL hybrid inversion implemented in the package prosail. \label{fig:hybrid}](FlowChart_HybridInversion_JOSS.tif){ width=75% }


### Application to Sentinel-2 image and comparison with SNAP

Once ML regression models are trained, `prosail` hybrid inversion can be applied
on Sentinel-2 imagery using the function `apply_prosail_inversion`.

The biophysical variables produced with SNAP are compared to those produced
with the `prosail` hybrid inversion.
The Sentinel-2 Level-2A product corresponding to tile *30SWJ* acquired on
*May 13th, 2021* was downloaded from the Copernicus Data Space Ecosystem
([CDSE](https://browser.dataspace.copernicus.eu)).


LAI, fcover and fAPAR computed from Sentinel-2 images with the hybrid inversion 
implemented in `prosail` are compared to those produced with SNAP in Figure 
\ref{fig:SNAP_prosail}. 

![Biophysical properties estimated with `prosail` vs Biophysical properties estimated with SNAP. \label{fig:SNAP_prosail}](compare_SNAP_prosail.png){ width=90% }

The two methods show good consistency, with Pearson correlation 
coefficient ranging between 0.99 and 1.00 for LAI, and superior to 0.97 for 
fCover and fAPAR.
Algorithmic differences remain between the results obtained from the biophysical 
processor of the SNAP toolbox and the inversion obtained with `prosail`, 
including: 
- the version of the PROSPECT model: SNAP uses a version anterior to PROSPECT-D,
- the ML algorithm. 

The differences in the retrieval of the vegetation biophysical properties suggest
differences in the soil properties accounted for in simulations with low LAI, 
despite our efforts to reproduce the workflow described in the ATBD [@weiss2020]. 
Additional tests performed over other sites (on both croplands and forests) 
showed similar performances. 

The availability of open source and fully parameterizable inversion procedures 
should contribute to improve reproducibility of currently available softwares.

# Conclusion

We introduce `prosail`, an R package dedicated to the canopy reflectance model 
PROSAIL. 
`prosail` is coupled with the R package `prospect` in order to allow 
seamless integration of future versions of the leaf model. 
It allows simulation of bi-directional reflectance factor and surface 
reflectance for any type of optical sensor based on their spectral response. 
The package also includes a collection of inversion procedures, including 
iterative optimization with and without prior information, and hybrid inversion. 

The estimation of vegetation biophysical properties with `prosail` hybrid 
inversion is consistent with estimations from the reference Sentinel toolbox. 

`prosail` hybrid inversion does not intend to be computationally as efficient as 
SNAP and may not meet requirements for large scale vegetation monitoring 
applications. 
`prosail` hybrid inversion offers a fully adjustable hybrid inversion framework, 
including an original parsimonious ML regression method, allowing 
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
We are grateful to Stéphane Jacquemoud and Frédéric Baret for the development 
of the initial version of the prospect model.

# AI Usage Disclosure

No generative AI tools were used in the development of this software, the 
writing of this manuscript, or the preparation of supporting materials.

# References
