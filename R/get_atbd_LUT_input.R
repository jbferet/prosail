#' This function provides a LUT corresponding to PROSAIL input parameters to
#' reproduce the distribution of parameters used in the ATBD
#'
#' @param nbSamples numeric. number of samples to generate
#' @param GeomAcq list. min and max values for tts, tto and psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#'
#' @return InputPROSAIL list. list of PROSAIL input parameters
#' @importFrom truncnorm rtruncnorm
#' @export

get_atbd_LUT_input <- function(nbSamples = 2000, GeomAcq = NULL,
                               Codist_LAI = TRUE){
  # define paremertization for truncated gaussians
  TgaussParms <- list()
  TgaussParms$min <- data.frame('lai' = 0, 'LIDFa' = 30, 'q' = 0.1, 'N' = 1.2,
                                'CHL' = 20, 'LMA' = 0.003, 'Cw_rel' = 0.6,
                                'BROWN' = 0.0, 'psoil' = 0)
  TgaussParms$max <- data.frame('lai' = 15, 'LIDFa' = 80, 'q' = 0.5, 'N' = 1.8,
                                'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.85,
                                'BROWN' = 2.0, 'psoil' = 1)
  TgaussParms$mean <- data.frame('lai' = 2, 'LIDFa' = 60, 'q' = 0.2, 'N' = 1.5,
                                 'CHL' = 45, 'LMA' = 0.005, 'Cw_rel' = 0.75,
                                 'BROWN' = 0.0, 'psoil' = 0.25)
  TgaussParms$sd <- data.frame('lai' = 3, 'LIDFa' = 30, 'q' = 0.5, 'N' = 0.3,
                               'CHL' = 30, 'LMA' = 0.005, 'Cw_rel' = 0.08,
                               'BROWN' = 0.30, 'psoil' = 0.6)

  # get distribution corresponding to gaussians
  InputPROSAIL <- list()
  for (parm in names(TgaussParms$min)){
    set.seed(42)
    InputPROSAIL[[parm]] <- truncnorm::rtruncnorm(n = nbSamples,
                                                  a = TgaussParms$min[[parm]],
                                                  b = TgaussParms$max[[parm]],
                                                  mean = TgaussParms$mean[[parm]],
                                                  sd = TgaussParms$sd[[parm]])
  }
  InputPROSAIL <- data.frame(InputPROSAIL)

  # define co-distribution with LAI
  if (Codist_LAI==TRUE){
    Codist_LAI <- list()
    Codist_LAI$Vmin0 <- data.frame('LIDFa' = 30, 'q' = 0.1, 'N' = 1.2,
                                   'CHL' = 20, 'LMA' = 0.003, 'Cw_rel' = 0.6,
                                   'BROWN' = 0.0, 'psoil' = 0)
    Codist_LAI$Vmax0 <- data.frame('LIDFa' = 80, 'q' = 0.5, 'N' = 1.8,
                                   'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.85,
                                   'BROWN' = 2.0, 'psoil' = 1)
    Codist_LAI$VminLAImax <- data.frame('LIDFa' = 55, 'q' = 0.1, 'N' = 1.3,
                                        'CHL' = 45, 'LMA' = 0.005,
                                        'Cw_rel' = 0.70,'BROWN' = 0.0,
                                        'psoil' = 0)
    Codist_LAI$VmaxLAImax <- data.frame('LIDFa' = 65, 'q' = 0.5, 'N' = 1.8,
                                        'CHL' = 90, 'LMA' = 0.011,
                                        'Cw_rel' = 0.80, 'BROWN' = 0.2,
                                        'psoil' = 0.4)
    for (parm in names(Codist_LAI$Vmin0)){
      Vstar <- get_codistributions(V = InputPROSAIL[[parm]],
                                   LAI = InputPROSAIL$lai,
                                   MaxLAI = TgaussParms$max$lai,
                                   Vmin0 = Codist_LAI$Vmin0[[parm]],
                                   Vmax0 = Codist_LAI$Vmax0[[parm]],
                                   VminLAImax = Codist_LAI$VminLAImax[[parm]],
                                   VmaxLAImax = Codist_LAI$VmaxLAImax[[parm]])
      InputPROSAIL[[parm]] <- Vstar
    }
  }

  # convert CW_rel into EWT
  InputPROSAIL$EWT <- ((InputPROSAIL$LMA)/(1-InputPROSAIL$Cw_rel)) -
    InputPROSAIL$LMA
  # set ANT, PROT and CBC to 0
  InputPROSAIL$ANT <- InputPROSAIL$PROT <- InputPROSAIL$CBC <- 0
  # set CAR to 0.25*CHL
  InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL
  # set geometry of acquisition
  if (is.null(GeomAcq))
    GeomAcq <- data.frame('min' = c('tto' = 0, 'tts' = 20, 'psi' = 0),
                          'max' = c('tto' = 10, 'tts' = 30, 'psi' = 360))
  if (inherits(GeomAcq, "list"))
    GeomAcq <- data.frame('min' = c('tto' = GeomAcq$min$tto,
                                    'tts' = GeomAcq$min$tts,
                                    'psi' = GeomAcq$min$psi),
                          'max' = c('tto' = GeomAcq$max$tto,
                                    'tts' = GeomAcq$max$tts,
                                    'psi' = GeomAcq$max$psi))
  InputPROSAIL$tts <- runif(n = nbSamples, min = GeomAcq['tts', 'min'],
                            max = GeomAcq['tts', 'max'])
  InputPROSAIL$tto <- runif(n = nbSamples, min = GeomAcq['tto', 'min'],
                            max = GeomAcq['tto', 'max'])
  InputPROSAIL$psi <- runif(n = nbSamples, min = GeomAcq['psi', 'min'],
                            max = GeomAcq['psi', 'max'])
  # default values
  InputPROSAIL$TypeLidf <- 2
  InputPROSAIL$alpha <- 40
  return(InputPROSAIL)
}
