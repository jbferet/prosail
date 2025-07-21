#' generates a LUT of brf based on a table of input variables for prosail model
#'
#' @param input_prosail list. PROSAIL input variables
#' @param spec_prospect list. Includes optical constants required for PROSPECT
#' @param spec_soil list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param spec_atm list. Includes direct and diffuse radiation for
#' clear conditions
#' @param band_names character. Name of the spectral bands of the sensor
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param brown_lop list. Defines optical properties for brown vegetation
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#'
#' @return brf numeric. matrix of brf corresponding to input_prosail
#' @importFrom progress progress_bar
#' @export

generate_lut_brf <- function(input_prosail, spec_prospect, spec_soil, spec_atm,
                             band_names = NULL, SAILversion = '4SAIL',
                             brown_lop = NULL){

  nb_samples <- length(input_prosail[[1]])
  brf <- list()
  split_nb <- round(nb_samples/10)
  pb <- progress_bar$new(
    format = "Generate LUT [:bar] :percent in :elapsed",
    total = 10, clear = FALSE, width= 100)
  for (i in seq_len(nb_samples)){
    if (i%%split_nb==0 & nb_samples>100) pb$tick()
    rsoil <- input_prosail[i,]$psoil*spec_soil$max_refl +
      (1-input_prosail[i,]$psoil)*spec_soil$min_refl
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      refl <- prosail(spec_sensor = spec_prospect,
                      input_prospect = input_prosail[i,],
                      type_lidf = input_prosail[i,]$type_lidf,
                      lidf_a = input_prosail[i,]$lidf_a,
                      lidf_b = input_prosail[i,]$lidf_b,
                      lai = input_prosail[i,]$lai,
                      hotspot = input_prosail[i,]$hotspot,
                      tts = input_prosail[i,]$tts,
                      tto = input_prosail[i,]$tto,
                      psi = input_prosail[i,]$psi, rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      refl <- prosail(spec_sensor = spec_prospect,
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
    # Computes brf based on outputs from PROSAIL and sun position
    brf[[i]] <- compute_brf(rdot = refl$rdot,
                            rsot = refl$rsot,
                            tts = input_prosail$tts[[i]],
                            spec_atm_sensor = spec_atm)
  }
  brf <- do.call(cbind,brf)
  row.names(brf) <- band_names
  return(brf)
}
