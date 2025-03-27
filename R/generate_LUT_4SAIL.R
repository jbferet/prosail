#' generates a LUT of 4SAIL outputs based on a table of input variables for PRO4SAIL model
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
#'
#' @return 4SAIL_LUT numeric. matrix of 4SAIL outputs (rdot, rsot, rsdt, rddt)
#' corresponding to input_prosail
#'
#' @importFrom progress progress_bar
#' @export

generate_LUT_4SAIL <- function(input_prosail, SpecPROSPECT, SpecSOIL, SpecATM,
                               band_names = NULL, SAILversion ='4SAIL',
                               brown_lop = NULL){

  nb_samples <- length(input_prosail[[1]])
  rdot <- rsot <- rsdt <- rddt <-  BRF <- list()
  Split <- round(nb_samples/10)
  pb <- progress_bar$new(
    format = "Generate LUT [:bar] :percent in :elapsed",
    total = 10, clear = FALSE, width= 100)
  for (i in seq_len(nb_samples)){
    if (i%%Split==0 & nb_samples>100) pb$tick()
    rsoil <- input_prosail[i,]$psoil*SpecSOIL$Dry_Soil +
      (1-input_prosail[i,]$psoil)*SpecSOIL$Wet_Soil
    # if 4SAIL
    if (SAILversion=='4SAIL'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          input_prospect = input_prosail[i,],
                          TypeLidf = input_prosail[i,]$TypeLidf,
                          LIDFa = input_prosail[i,]$LIDFa,
                          LIDFb = input_prosail[i,]$LIDFb,
                          lai = input_prosail[i,]$lai,
                          q = input_prosail[i,]$q,
                          tts = input_prosail[i,]$tts,
                          tto = input_prosail[i,]$tto,
                          psi = input_prosail[i,]$psi,
                          rsoil = rsoil)
    } else if (SAILversion=='4SAIL2'){
      RefSAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                          input_prospect = input_prosail[i,],
                          TypeLidf = input_prosail[i,]$TypeLidf,
                          LIDFa = input_prosail[i,]$LIDFa,
                          LIDFb = input_prosail[i,]$LIDFb,
                          lai = input_prosail[i,]$lai,
                          q = input_prosail[i,]$q,
                          tts = input_prosail[i,]$tts,
                          tto = input_prosail[i,]$tto,
                          psi = input_prosail[i,]$psi, rsoil = rsoil,
                          SAILversion = '4SAIL2',
                          fraction_brown = input_prosail[i,]$fraction_brown,
                          diss = input_prosail[i,]$diss,
                          Cv = input_prosail[i,]$Cv,
                          Zeta = input_prosail[i,]$Zeta, brown_lop = brown_lop)
    }
    rdot[[i]] <- RefSAIL$rdot
    rsot[[i]] <- RefSAIL$rsot
    rsdt[[i]] <- RefSAIL$rsdt
    rddt[[i]] <- RefSAIL$rddt
    # Computes BRF based on outputs from PROSAIL and sun position
    BRF[[i]] <- compute_BRF(rdot = RefSAIL$rdot, rsot = RefSAIL$rsot,
                            tts = input_prosail$tts[[i]],
                            SpecATM_Sensor = SpecATM)
  }
  BRF <- do.call(cbind,BRF)
  rdot <- do.call(cbind,rdot)
  rsot <- do.call(cbind,rsot)
  rsdt <- do.call(cbind,rsdt)
  rddt <- do.call(cbind,rddt)
  row.names(BRF) <- row.names(rdot) <- row.names(rsot) <- band_names
  row.names(rsdt) <- row.names(rddt) <- band_names
  return(list('BRF' = BRF, 'rdot' = rdot, 'rsot' = rsot,
              'rsdt' = rsdt, 'rddt' = rddt))
}
