#' This function generates a LUT of PROSAIL outputs, including BRF, fAPAR,
#' fCover and albedo, based on a table of input variables for PRO4SAIL model
#'
#' @param input_prosail list. PROSAIL input variables
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
#' @param band_names character. Name of the spectral bands of the sensor
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param brown_lop list. Defines optical properties for brown vegetation
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#' @param progress boolean. set TRUE for progress bar during production of LUT
#'
#' @return LUT numeric. list of BRF, fCover, fAPAR and albedo corresponding
#' to input_prosail
#' @importFrom progress progress_bar
#' @export

generate_LUT_PROSAIL <- function(input_prosail, SpecPROSPECT, SpecSOIL, SpecATM,
                                 band_names = NULL, SAILversion = '4SAIL',
                                 brown_lop = NULL, progress = TRUE){

  nb_samples <- length(input_prosail[[1]])
  BRF <- list()
  fCover <- fAPAR <- albedo <- c()
  Split <- round(nb_samples/10)
  if (progress)
    pb <- progress_bar$new(
      format = "Generate LUT [:bar] :percent in :elapsed",
      total = 10, clear = FALSE, width= 100)
  for (i in seq_len(nb_samples)){
    if (progress & i %% Split == 0 & nb_samples>100) pb$tick()
    rsoil <- input_prosail[i,]$psoil*SpecSOIL$Dry_Soil +
      (1-input_prosail[i,]$psoil)*SpecSOIL$Wet_Soil
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          input_prospect = input_prosail[i,],
                          TypeLidf = input_prosail$TypeLidf[i],
                          LIDFa = input_prosail[i,]$LIDFa,
                          LIDFb = input_prosail[i,]$LIDFb,
                          lai = input_prosail[i,]$lai,
                          q = input_prosail[i,]$q,
                          tts = input_prosail[i,]$tts,
                          tto = input_prosail[i,]$tto,
                          psi = input_prosail[i,]$psi, rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          input_prospect = input_prosail[i,],
                          TypeLidf = input_prosail[i,]$TypeLidf,
                          LIDFa = input_prosail[i,]$LIDFa,
                          LIDFb = input_prosail[i,]$LIDFb,
                          lai = input_prosail[i,]$lai, q = input_prosail[i,]$q,
                          tts = input_prosail[i,]$tts,
                          tto = input_prosail[i,]$tto,
                          psi = input_prosail[i,]$psi, rsoil = rsoil,
                          SAILversion = '4SAIL2',
                          fraction_brown = input_prosail[i,]$fraction_brown,
                          diss = input_prosail[i,]$diss,
                          Cv = input_prosail[i,]$Cv,
                          Zeta = input_prosail[i,]$Zeta,
                          brown_lop = brown_lop)
    }
    # Computes BRF based on outputs from PROSAIL and sun position
    BRF[[i]] <- compute_BRF(rdot = RefSAIL$rdot,
                            rsot = RefSAIL$rsot,
                            tts = input_prosail$tts[[i]],
                            SpecATM_Sensor = SpecATM)
    fCover[i] <- RefSAIL$fCover
    fAPAR[i] <- compute_fAPAR(abs_dir = RefSAIL$abs_dir,
                              abs_hem = RefSAIL$abs_hem,
                              tts = input_prosail$tts[[i]],
                              SpecATM_Sensor = SpecATM)
    albedo[i] <- compute_albedo(rsdstar = RefSAIL$rsdstar,
                                rddstar = RefSAIL$rddstar,
                                tts = input_prosail$tts[[i]],
                                SpecATM_Sensor = SpecATM)
  }
  BRF <- do.call(cbind,BRF)
  row.names(BRF) <- band_names
  return(list('BRF' = BRF, 'fCover' = fCover,
              'fAPAR' = fAPAR, 'albedo' = albedo))
}
