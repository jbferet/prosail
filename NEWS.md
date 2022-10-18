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