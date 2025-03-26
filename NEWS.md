# prosail v2.5.10

## addition
- added SRF for Sentinel-2C and updated SRF for 2A and 2B
- added test-PRO4SAIL_iterative_opt

## change
- create directories of not created when training and applying inversion

# prosail v2.5.9

## change
- added skyl as input parameter for Compute_BRF

# prosail v2.5.8

## change
- added option for progress bar in Generate_LUT_PROSAIL

# prosail v2.5.7

## change
- added '.tiff' extension when needed when setting 'bigRaster = T' in Apply_prosail_inversion

# prosail v2.5.6

## fix
- correction on PROSAIL_Hybrid_Apply: seq_len(nbEnsemble) instead of seq_len(length(nbEnsemble))


# prosail v2.5.5

## addition
- possibility to select file format for output BP rasters. default = 'GTiff' 

# prosail v2.5.4

## addition
- finalize introduction of 4SAIL2 in PROSAIL inversion: introduce BrownLOP as 
input variable for functions 'train_prosail_inversion'. 'Generate_LUT_PROSAIL' 
and 'Generate_LUT_BRF'

# prosail v2.5.3

## fix
- PROSAIL_Hybrid_Train: convert BRF_LUT into data.frame to apply dplyr::slice

## change
- update tutorial #2: 'wvl <- SpecPROSPECT_FullRange$lambda'

# prosail v2.5.2

## fix
- use dplyr slice to select rows from dataframe when performing hybrid inversion 
(allows using 1 variable only with liquidSVM)

# prosail v2.5.1

## addition
- added option for progressbar

# prosail v2.5.0

## addition
- supports bigRaster for the application of PROSAIL hybrid inversion on raster data
- implement additional ML algorithms for inversion. Currently suboptimal

# prosail v2.4.1

## fix
- fixed a bug occurring when running 4SAIL2 with brown vegetation and Input_PROSPECT undefined

# prosail v2.4.0

## fix
- fixed the tutorial related to hybrid inversion applied on Sentinel-2 image: B8A was not listed in spectral bands for the raster

# prosail v2.3.3

## change
- move liquidSVM from Imports to Suggests
- remove unnecessary packages from Imports
- remove figure option for hybrid inversion training

# prosail v2.3.2

## fix
- convert InRefl in dataframe to make sure applySensorCharacteristics works even if 1 sample provided.
- fix Invert_PROSAIL

## changes
- simplify GetRadiometry 
- update read_ENVI_header to account for .hdr and .HDR

# prosail v2.3.1

## additions
- added automated unit tests
- include plain text versions for SpecPROSPECT, SpecATM and SpecSOIL 

## changes
- changed default values

# prosail v2.3.0

## Changes
- update to fix bugs induced by upgrade of package prospect (>=1.6.0)
- rewrite PRO4SAIL function to simplify it

# prosail v2.2.3

## Changes
- added function Generate_LUT_4SAIL to produce LUT corresponding to all 4SAIL outputs (rdot, rsot, rsdt, rddt) in addition to BRF

# prosail v2.2.2

## Changes
- change format of GeomAcq in get_atbd_LUT_input: using data.frame instead of list

# prosail v2.2.1

## Fix
- correct Apply_prosail_inversion so that inversion can be performed on images with file extension

# prosail v2.2.0

## Changes
- add fAPAR, albedo and fCover as possible variables to estimate from train_prosail_inversion

# prosail v2.1.0

## Changes
- add function Compute_fAPAR and Compute_albedo to compute fAPAR and albedo

# prosail v2.0.2

## Changes
- eliminate LMA when generating PROSAIL-PRO parameter distribution with get_distribution_input_prosail in order to avoid warnings when producing reflectance
- function PROSAIL() : sets LMA to 0 when LMA = NULL

# prosail v2.0.1

## Changes
- add function apply_noise_AddMult to apply additive and multiplicative noise to reflectance LUT

# prosail v2.0.0

## Changes
- added possibility to use input PROSAIL distribution and noise model from ATBD when performing hybrid inversion
- added possibility to use user defined Input PROSAIL and user defined BRF LUT
- updated vignettes
- removed 4SAIL2 from option for hybrid inversion, unless users provide their own input parameter distribution and BRF LUT

# prosail v1.3.2

## Fix
- fixed bug occuring when atbd == NULL


# prosail v1.3.1:

## Fix
- fixed train_prosail_inversion : systematically add name for spectral bands in BRF_LUT according to SRF when providing


# prosail v1.3.0

## Changes
- implemented configuration described in S2 toolbox ATBD for biophysical processor (http://step.esa.int/docs/extra/ATBD_S2ToolBox_V2.1.pdf)
 - truncated gaussian distribution
 - co-distribution with LAI
 - additive and multiplicative noise
- added possibility to directly provide user-defined InputPROSAIL data.frame to complement ATBD and initial minmax range and distributions
- added spectral response for Landsat-9 & Pleiades, added central wavelength and sensor name in SRF variable
- returning path for biophysical map products in Apply_prosail_inversion
- removed dependency to rgdal
- function added: apply_noise_atbd, get_atbd_LUT_input, get_codistributions, get_default_LUT_input

# prosail v1.2.4

## Changes
- added raster file path as output of function Apply_prosail_inversion. corresponds to maps of biophysical properties of interest (mean value and corresponding standard deviation)

# prosail v1.2.3

## Changes
- added function OptimalSI to compute the correlation between a set of vegetation properties and all combinations of spectral bands corresponding to a given type of spectral index

## Fix
- transpose sensor spectral response function in applySensorCharacteristics if not properly oriented
- Apply_Noise_LUT: add possibility to add absolute noise

# prosail v1.2.2

## Fix
- corrected function to compute spectral indices: SpectralIndices instead of spectralindices
- set default value for LMA to 0

# prosail v1.2.1

## Changes
- implemented caret for function Apply_prosail_inversion used for raster processing
- added verbose parameter to control printed message when hyperparameter adjustment performed during training of liquidSVM

# prosail v1.2.0

## Changes
- added SVM from caret as alternative to liquidSVM

# prosail v1.1.1

## Changes
- added CR_RE in spectral indices

# prosail v1.0.0

## Changes
- application of hybrid inversion on rasters
- updated spectral indices computed from reflectance matrix
- added function to get S2 geometry from MTD_TL.xml files

# prosail v0.9.0

## Changes
- added hybrid inversion based on liquidSVM
- created Lib_PROSAIL_LUT
- updated spectral indices computed from reflectance matrix
- update email addresses
- modified default diss value (4SAIL2) to 0
- modified default fraction_brown value (4SAIL2) to 0 and allow only one leaf optics to be defined when fraction_brown = 0

## Fixes
- refined import to avoid warnings

# prosail v0.1.0

## Fixes
- fix SOILSPECT input parameters

## Changes
- introduces 4SAIL2 in addition to 4SAIL
- directly calls the PROSPECT package instead of including the source code
- added Lib_SpectralIndices.R
- SOILSPECT added in Lib_SOILSPECT.R

# prosail v0.0.1

## Fixes
- fix SOILSPECT input parameters
- update PROSPECT for handling of error during inversion
- modified instructions fo install when prosail not public
- fix dladgen in Lib_PROSAIL.R
- fix warning of dladgen in Lib_PROSAIL.R
- code cleaning in Lib_PROSAIL.R

## Changes
- introduced iterative inversion
- added prior information to inversion
- started working on vignettes
- Upgrade version in DESCRIPTION
- Include Venus and Sentinel-2 sensors
