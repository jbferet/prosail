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

PROSAIL is a vegetation radiative transfer model (RTM) coupling the 
*leaf optical PROperties SPECtra* (PROSPECT) model, which simulates leaf optical 
properties [@jacquemoud1990], with the *Scattering by Arbitrarily Inclined Leaves* 
model (SAIL) which simulates canopy bidirectional reflectance 
[@verhoef1984; @jacquemoud2009; @berger2018] from the visible to shortwave 
infrared domains. 
PROSAIL is a valuable tool for agriculture, ecology, and climate research [@poulter2023]. 
Inversion procedures allow retrieval of vegetation biophysical properties, 
including 

- leaf area index (LAI)

- leaf chemistry

- the fraction of absorbed photosynthetically active radiation (fAPAR)

- fractional vegetation cover (fCover)

The retrieval accuracy of these properties depends on the spectral 
characteristics of the sensor: solving an ill-posed problem requires sufficient 
spectral information. 

We introduce `prosail`, an R package which provides multiple versions of 
the model PROSAIL, coupled `prospect` [@feret2024]. 
`prosail` simulates surface reflectance in forward mode for various optical 
sensors, based on their spectral response function (SRF).
It includes functions for model inversion using iterative optimization, and a 
hybrid approach, combining physical modeling and machine learning (ML) regression. 

# Statement of need

Vegetation traits mapping and monitoring is crucial to better understand 
ecosystem functions, carbon, water and energy budgets. 
Vegetation RTMs describe the interactions between light and vegetation
based on their biophysical and chemical properties, through absorption and 
scattering mechanisms. 

The SAIL model is an example of four-stream representations of the radiative 
transfer equation, including two direct fluxes (incident solar flux and radiance 
in the viewing direction) and two diffuse fluxes (upward and downward 
hemispherical flux).

More complex models include the *Soil Canopy Observation of Photosynthesis and 
Energy fluxes* (SCOPE) [@yang2021] model, which features layered description of 
the vegetation, and the *Discrete Anisotropic Radiative Transfer* DART model 
[@gastellu2015], which explicitly describes 3D canopy structure.

Hybrid strategies generating a training dataset for ML algorithms with RTM 
simulations show multiple advantages [@verrelst2015]: simulations complement or 
replace in situ sampling for robust estimation, while ML ensures computation 
efficiency compared to iterative optimization.

The limited number of input parameters and the computational efficiency of SAIL
is interesting for operational applications, despite the main assumption of 
vegetation as a homogeneous turbid medium: leaves are randomly distributed 
within the canopy, which is not accurate for row crops and heterogeneous canopies.

# State of the field

Various softwares allow hybrid inversion with PROSAIL simulations. 
The *Sentinel Toolbox Application Platform* (SNAP) includes the 
`Biophysical Processor` module [@weiss2020], combining PROSAIL simulations with 
an artificial neural network regression model. 
The sampling design of the training set is described in the 
*Algorithmic Theoretical Background Document* (ATBD) [@weiss2020]).
Alternative distributions provide more interactive parameterization of the 
inversion strategy, e.g. the *Automated Radiative Transfer Models Operator*
(ARTMO) Matlab toolbox [@rivera2014].

`prosail` includes *4SAIL* [@verhoef2007], and *4SAIL2* [@verhoefbach2007], 
a two layers version of *4SAIL*. 

`prosail` does not intend to provide the same computational efficiency as SNAP. 
It does not provide a collection of models and methods as comprehensive as those 
provided with the ARTMO box neither. 
`prosail` provides a flexible open source framework to experiment with hybrid 
inversion procedures, using simple yet fast and efficient training stage, 
allowing experimenting on training data sampling design, introduction 
of noise, feature selection over optical sensors. 
It is suitable for research and development, to identify potential and 
limitations of inversion strategies.
Resulting regression models can be applied on remote sensing data, but it may 
not be the appropriate software for scaling up operational vegetation monitoring 
applications.

Alternative PROSAIL implementations are also available in 
[python](https://github.com/jgomezdans/prosail), 
[Julia](https://github.com/RemoteSensingTools/CanopyOptics.jl) and 
[R](https://github.com/ashiklom/rrtm).
Other options are available on 
[this webpage](http://teledetection.ipgp.jussieu.fr/prosail/).


# Software design

`prosail` uses `prospect` for the simulation of leaf optical properties. 
`prosail` provides user-friendly and modular functions to simulate vegetation 
canopy reflectance, and to predict biophysical properties using inversion. 

Hybrid inversion is a multi-step procedure requiring simulation of sensor 
reflectance to train a ML regression algorithm, then applicable to any data 
source collected from airborne and spaceborne sensors. 

- The generation of realistic reflectance simulations requires definition of 
input parameter sampling strategy, accounting for sensor SRF and uncertainties. 
`prosail` provides functions to produce this training dataset as a unique step 
following standard parameterization, or step by step, in order for users to 
understand the logic, or adjust each step to their needs.

- A limited number of ML algorithms is currently provided.
Additional algorithms will be implemented in future versions.

- Raster processing is handled with the `terra` [@terra] package, for 
compatibility with most raster data formats. 

- `prosail` can be used in combination with 
[`preprocS2`](https://jbferet.gitlab.io/preprocs2/), a package dedicated to 
downloading and preprocessing of Earth observation data following the 
[STAC](https://stacspec.org/en) specification, in order to produce a fully 
automated image access and processing workflow.

# Research impact statement

`prosail` has been used in multiple research publications since its early 
developments [@hauser2021], [@ferreira2026], [@feret2026], [@kattenborn_temporal_2024].
It is currently used in multiple research projects, and is actively maintained. 


# Overview

## PROSAIL simulation in forward mode

PROSAIL requires information intrinsic to vegetation : 

- leaf optical properties, either simulated with `prospect`, or measured
experimentally.

- LAI

- Leaf inclination distribution, defined from different options of 
*Leaf Inclination Distribution Functions* (LIDF). 

- foliage hot spot parameter [@kuusk1991; @breon2002; @verhoefbach2007]

PROSAIL also requires information extrinsic to vegetation : 

- sun-observer geometry: sun and observer zenith angles, relative azimuth angle

- soil reflectance

SAIL produces top-of-canopy bi-hemispherical reflectance factor, 
directional-hemispherical reflectance factor for solar incident flux, 
hemispherical-directional reflectance factor in viewing direction and 
bi-directional reflectance factor.


## Simulating surface reflectance acquired by a sensor 

Simulating sensor surface reflectance requires a SRF. 
A selection of SRF corresponding to multiple sensors is already implemented in 
`prosail`, including Sentinel-2, Landsat and MODIS.

User-defined sensors require defining central wavelength and full width at half 
maximum (fwhm) for each band, assuming gaussian response. 
Exact SRF from external files can also be provided. 


## PROSAIL hybrid inversion

### Simulating surface reflectance to prepare for regression model training 

Preparing for training regression models needs defining: 

- a sampling strategy (range, distribution, co-distributions for input 
parameters)

- a level of noise applicable to the simulated surface reflectance

- spectral bands used for the training of the regression model

The default parameterization of `prosail` hybrid inversion relies on the 
production of a surface reflectance LUT compliant with the specifications 
defined in the ATBD for the Biophysical Processor of the Sentinel toolbox [@weiss2020]. 
Each step of the simulation of a training LUT can be defined by user.


### Training a ML regression algorithm

The ML strategy implemented in `prosail` is a parsimonious ensemble method based 
on a bootstrap aggregating (bagging) prediction of biophysical properties from 
support vector regression (SVR) models.
Each SVR model is trained with a limited number of samples, ensuring fast 
training stage. 
The predicted value then corresponds to the mean prediction from a set of 
SVR models.
The standard deviation is derived from this ensemble of predictors.
However, standard deviation may not provide accurate model uncertainty 
quantification [@palmer2022].

The default parameterization uses an ensemble of 10 SVR models trained with 200 
samples each. 
The training then requires 2000 samples, to be compared to the 41472 simulations
used for training models implemented in SNAP. 
This parsimonious ensemble method ensures fast training stage and retrieval 
performances similar to SNAP. 

### Performing full hybrid inversion

Figure \ref{fig:hybrid} summarizes the workflow applied to perform hybrid 
inversion. 
The function `train_prosail_inversion` combines aforementioned steps to 
simulate surface reflectance and train SVR models. 
Obtained regression models are then used to estimate vegetation biophysical 
properties from data tables or raster data. 

![Workflow of PROSAIL hybrid inversion implemented in the package prosail. \label{fig:hybrid}](FlowChart_HybridInversion_JOSS.tif){ width=75% }


### Application to Sentinel-2 image and comparison with SNAP

Once ML regression models are trained, `prosail` hybrid inversion can be applied
on Sentinel-2 imagery using the function `apply_prosail_inversion`.

Biophysical variables produced with SNAP are compared to those produced
with the `prosail` hybrid inversion.
The Sentinel-2 Level-2A product corresponding to tile *30SWJ* acquired on
*May 13th, 2021* was downloaded from the Copernicus Data Space Ecosystem
([CDSE](https://browser.dataspace.copernicus.eu)).


LAI, fcover and fAPAR computed from Sentinel-2 images with `prosail` hybrid 
inversion are compared to those produced with SNAP in Figure \ref{fig:SNAP_prosail}. 

![Biophysical properties estimated with `prosail` vs Biophysical properties estimated with SNAP. \label{fig:SNAP_prosail}](compare_SNAP_prosail.png){ width=90% }

The two methods show good consistency, with Pearson correlation 
coefficient > 0.99 for LAI, and > 0.97 for fCover and fAPAR.
Algorithmic differences remain between the two implementations, including: 
- the version of the PROSPECT model: SNAP uses a version anterior to PROSPECT-D,
- the ML algorithm. 

Differences in soil properties accounted for in simulations with low LAI may 
also contribute, despite efforts to reproduce the workflow described in the ATBD. 
Additional tests performed over croplands and forests showed similar performances. 

The availability of open source and fully parameterizable inversion procedures 
should contribute to improve reproducibility of currently available softwares.

# Conclusion

We introduce `prosail`, an R package dedicated to the canopy reflectance model 
PROSAIL. 
`prosail` simulates canopy reflectance for optical sensors based on their 
spectral response. 
The package includes inversion based on iterative optimization and ML/RTM 
inversion. 

The estimation of vegetation biophysical properties with `prosail` hybrid 
inversion is consistent with estimations from SNAP. 

`prosail` hybrid inversion does not intend to be computationally as efficient as 
SNAP and may not meet requirements for large scale vegetation monitoring 
applications. 
`prosail` hybrid inversion offers a fully adjustable hybrid inversion framework, 
including an original parsimonious ML regression method, allowing fast and 
efficient training. 
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
