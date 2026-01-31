#' This function generates a LUT of PROSAIL outputs, including surface reflectance, fapar,
#' fcover and albedo, based on a table of input variables for prosail model
#'
#' @param input_prosail list. PROSAIL input variables
#' @param spec_prospect list. Includes optical constants required for PROSPECT
#' @param spec_soil list. Includes either a set of OSSL library reflectance
#' spectra, reflectance spectra from S2 ATBD v2, or minimum reflectance and maximum reflectance
#' @param spec_atm list. direct and diffuse radiation for clear conditions
#' @param band_names character. Name of the spectral bands of the sensor
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param brown_lop list. Defines optical properties for brown vegetation
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#' @param progress boolean. set TRUE for progress bar during production of LUT
#'
#' @return LUT numeric. list of surf_refl, fcover, fapar and albedo corresponding
#' to input_prosail
#' @importFrom progress progress_bar
#' @export

generate_lut_prosail <- function(input_prosail, spec_prospect, spec_soil,
                                 spec_atm, band_names = NULL,
                                 SAILversion = '4SAIL', brown_lop = NULL,
                                 progress = TRUE){

  nb_samples <- length(input_prosail[[1]])
  surf_refl <- list()
  fcover <- fapar <- albedo <- c()
  split_nb <- round(nb_samples/10)
  if (progress)
    pb <- progress_bar$new(
      format = "Generate LUT [:bar] :percent in :elapsed",
      total = 10, clear = FALSE, width= 100)
  for (i in seq_len(nb_samples)){

    if (progress & i %% split_nb == 0 & nb_samples>100)
      pb$tick()
    if (!is.null(input_prosail[i,]$soil_brightness) &
        !is.null(input_prosail[i,]$soil_ID)){
      rsoil <- input_prosail[i,]$soil_brightness*
        spec_soil[[input_prosail[i,]$soil_ID+1]]
    } else if (is.null(input_prosail[i,]$soil_brightness) &
               ! is.null(input_prosail[i,]$psoil)){
      rsoil <- input_prosail[i,]$psoil*spec_soil$max_refl +
        (1-input_prosail[i,]$psoil)*spec_soil$min_refl
    }
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      refl_sail <- prosail(spec_sensor = spec_prospect,
                           input_prospect = input_prosail[i,],
                           type_lidf = input_prosail$type_lidf[i],
                           lidf_a = input_prosail[i,]$lidf_a,
                           lidf_b = input_prosail[i,]$lidf_b,
                           lai = input_prosail[i,]$lai,
                           hotspot = input_prosail[i,]$hotspot,
                           tts = input_prosail[i,]$tts,
                           tto = input_prosail[i,]$tto,
                           psi = input_prosail[i,]$psi, rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      refl_sail <- prosail(spec_sensor = spec_prospect,
                           input_prospect = input_prosail[i,],
                           type_lidf = input_prosail[i,]$type_lidf,
                           lidf_a = input_prosail[i,]$lidf_a,
                           lidf_b = input_prosail[i,]$lidf_b,
                           lai = input_prosail[i,]$lai,
                           hotspot = input_prosail[i,]$hotspot,
                           tts = input_prosail[i,]$tts,
                           tto = input_prosail[i,]$tto,
                           psi = input_prosail[i,]$psi, rsoil = rsoil,
                           SAILversion = '4SAIL2',
                           fraction_brown = input_prosail[i,]$fraction_brown,
                           diss = input_prosail[i,]$diss,
                           cv = input_prosail[i,]$cv,
                           zeta = input_prosail[i,]$zeta,
                           brown_lop = brown_lop)
    }
    # Computes surf_refl based on outputs from PROSAIL and sun position
    surf_refl[[i]] <- compute_surf_refl(rdot = refl_sail$rdot,
                                        rsot = refl_sail$rsot,
                                        tts = input_prosail$tts[[i]],
                                        spec_atm_sensor = spec_atm)
    fcover[i] <- refl_sail$fcover
    fapar[i] <- compute_fapar(abs_dir = refl_sail$abs_dir,
                              abs_hem = refl_sail$abs_hem,
                              tts = input_prosail$tts[[i]],
                              spec_atm_sensor = spec_atm)
    albedo[i] <- compute_albedo(rsdstar = refl_sail$rsdstar,
                                rddstar = refl_sail$rddstar,
                                tts = input_prosail$tts[[i]],
                                spec_atm_sensor = spec_atm)
  }
  surf_refl <- do.call(cbind,surf_refl)
  row.names(surf_refl) <- band_names
  return(list('surf_refl' = surf_refl, 'fcover' = fcover,
              'fapar' = fapar, 'albedo' = albedo))
}
