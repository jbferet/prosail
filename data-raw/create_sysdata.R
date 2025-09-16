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
SpecPROSPECT_FullRange <- read.table(file = 'data-raw/dataSpec_PRO.txt',
                                     header = TRUE,
                                     sep = '\t')

#' calctav_90 & calctav_40: transmissivity of a dielectric plane surface,
#' averaged over all directions of incidence and over all polarizations for
#' solid angle of 90 and 40 degrees
#'
calctav_90 <- prospect::calctav(90, nr = SpecPROSPECT_FullRange$nrefrac)
calctav_40 <- prospect::calctav(40, nr = SpecPROSPECT_FullRange$nrefrac)
SpecPROSPECT_FullRange$calctav_90 <- calctav_90
SpecPROSPECT_FullRange$calctav_40 <- calctav_40

# read dataspec file and save it in a dataframe
dataSoil_Atm <- read.table(file = 'data-raw/dataSpec_SOIL_ATM.txt',
                           header = TRUE,
                           sep = '\t')

SpecSOIL <- data.frame('lambda' = dataSoil_Atm$lambda,
                       'Dry_Soil' = dataSoil_Atm$Dry_Soil,
                       'Wet_Soil' = dataSoil_Atm$Wet_Soil)
SpecATM <- data.frame('lambda' = dataSoil_Atm$lambda,
                      'Direct_Light' = dataSoil_Atm$Direct_Light,
                      'Diffuse_Light' = dataSoil_Atm$Diffuse_Light)

usethis::use_data(SpecPROSPECT_FullRange, SpecSOIL, SpecATM,
                  internal = FALSE,
                  overwrite = TRUE)


spec_soil <- data.frame('lambda' = dataSoil_Atm$lambda,
                        'max_refl' = dataSoil_Atm$Dry_Soil,
                        'min_refl' = dataSoil_Atm$Wet_Soil)
spec_atm <- data.frame('lambda' = dataSoil_Atm$lambda,
                       'direct_light' = dataSoil_Atm$Direct_Light,
                       'diffuse_light' = dataSoil_Atm$Diffuse_Light)

spec_prospect_fullrange <- SpecPROSPECT_FullRange
usethis::use_data(spec_prospect_fullrange, spec_soil, spec_atm,
                  internal = FALSE,
                  overwrite = TRUE)


# read soil reflectance from Open Soil Spectral Library (OSSL)
# source #1: Sentinel-2 ATBD v3, 2023 https://hal.inrae.fr/hal-04884986v1/document
# source #2: Shepherd et al., 2022 https://doi.org/10.1016/j.soisec.2022.100061
spec_soil_ossl <- read.table(file = 'data-raw/dataSpec_SOIL_OSSL_ATBD_v3.txt',
                             header = TRUE,
                             sep = '\t')

usethis::use_data(spec_soil_ossl,
                  internal = FALSE,
                  overwrite = TRUE)
