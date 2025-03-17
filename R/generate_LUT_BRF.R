#' This function generates a LUT of BRF based on a table of input variables for
#' PRO4SAIL model
#'
#' @param InputPROSAIL list. PROSAIL input variables
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param SpecATM list. Includes direct and diffuse radiation for
#' clear conditions
#' @param BandNames character. Name of the spectral bands of the sensor
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param BrownLOP list. Defines optical properties for brown vegetation
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#'
#' @return BRF_LUT numeric. matrix of BRF corresponding to InputPROSAIL
#' @importFrom progress progress_bar
#' @export

generate_LUT_BRF <- function(InputPROSAIL, SpecPROSPECT, SpecSOIL, SpecATM,
                             BandNames = NULL, SAILversion='4SAIL',
                             BrownLOP = NULL){

  nbSamples <- length(InputPROSAIL[[1]])
  BRF <- list()
  Split <- round(nbSamples/10)
  pb <- progress_bar$new(
    format = "Generate LUT [:bar] :percent in :elapsed",
    total = 10, clear = FALSE, width= 100)
  for (i in seq_len(nbSamples)){
    if (i%%Split==0 & nbSamples>100) pb$tick()
    rsoil <- InputPROSAIL[i,]$psoil*SpecSOIL$Dry_Soil +
      (1-InputPROSAIL[i,]$psoil)*SpecSOIL$Wet_Soil
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          Input_PROSPECT = InputPROSAIL[i,],
                          TypeLidf = InputPROSAIL[i,]$TypeLidf,
                          LIDFa = InputPROSAIL[i,]$LIDFa,
                          LIDFb = InputPROSAIL[i,]$LIDFb,
                          lai = InputPROSAIL[i,]$lai,
                          q = InputPROSAIL[i,]$q,
                          tts = InputPROSAIL[i,]$tts,
                          tto = InputPROSAIL[i,]$tto,
                          psi = InputPROSAIL[i,]$psi, rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          Input_PROSPECT = InputPROSAIL[i,],
                          TypeLidf = InputPROSAIL[i,]$TypeLidf,
                          LIDFa = InputPROSAIL[i,]$LIDFa,
                          LIDFb = InputPROSAIL[i,]$LIDFb,
                          lai = InputPROSAIL[i,]$lai,
                          q = InputPROSAIL[i,]$q,
                          tts = InputPROSAIL[i,]$tts,
                          tto = InputPROSAIL[i,]$tto,
                          psi = InputPROSAIL[i,]$psi, rsoil = rsoil,
                          SAILversion = '4SAIL2',
                          fraction_brown = InputPROSAIL[i,]$fraction_brown,
                          diss = InputPROSAIL[i,]$diss,
                          Cv = InputPROSAIL[i,]$Cv,
                          Zeta = InputPROSAIL[i,]$Zeta,
                          BrownLOP = BrownLOP)
    }
    # Computes BRF based on outputs from PROSAIL and sun position
    BRF[[i]] <- compute_BRF(rdot = RefSAIL$rdot,
                            rsot = RefSAIL$rsot,
                            tts = InputPROSAIL$tts[[i]],
                            SpecATM_Sensor = SpecATM)
  }
  BRF <- do.call(cbind,BRF)
  row.names(BRF) <- BandNames
  return(BRF)
}
