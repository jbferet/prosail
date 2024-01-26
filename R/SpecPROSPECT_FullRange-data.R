#' SpecPROSPECT_FullRange: optical constants defined for PROSPECT
#'
#' Corresponds to spectral bands, refractive index and specific absorption
#' coefficient for each chemica constituent, defined over the spectral domain
#' ranging from 400 nm to 2500 nm
#'
#' The specific absorption coefficients were calibrated using experimental data.
#' The details of the calibration should be found in the following publications:
#'  http://dx.doi.org/10.1016/j.rse.2017.03.004
#'  https://doi.org/10.1016/j.rse.2020.112173
#'
#' @format A dataframe including 2101 rows and 10 columns
#'  - *lambda*: wavelength, between 400 nm and 2500 nm
#'  - *nrefrac*: refractive index
#'  - *SAC_CHL*: specific absorption coefficient for chlorophylls
#'  - *SAC_CAR*: specific absorption coefficient for carotenoids
#'  - *SAC_ANT*: specific absorption coefficient for anthocyanins
#'  - *SAC_BROWN*: specific absorption coefficient for brown pigments
#'  - *SAC_EWT*: specific absorption coefficient for water content (Equivalent Water Thickness)
#'  - *SAC_LMA*: specific absorption coefficient for dry matter (Leaf Mass per Area)
#'  - *SAC_PROT*: specific absorption coefficient for proteins
#'  - *SAC_CBC*: specific absorption coefficient for carbon-base constituents
#'  - *calctav_90*: transmissivity of a dielectric plane surface, averaged over
#'  all directions of incidence & all polarizations for solid angle of 90 deg
#'  - *calctav_40*: transmissivity of a dielectric plane surface, averaged over
#'  all directions of incidence & all polarizations for solid angle of 40 deg
"SpecPROSPECT_FullRange"
