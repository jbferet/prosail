#' This function provides a LUT corresponding to PROSAIL input parameters to
#' reproduce the distribution of parameters used in the ATBD
#'
#' @param nb_samples numeric. number of samples to generate
#' @param GeomAcq list. min and max values for tts, tto and psi
#' @param codistribution_lai boolean. set TYRUE if codistribution with LAI accounted for
#'
#' @return input_prosail list. list of PROSAIL input parameters
#' @importFrom truncnorm rtruncnorm
#' @export

get_atbd_LUT_input <- function(nb_samples = 2000, GeomAcq = NULL,
                               codistribution_lai = TRUE){
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
  input_prosail <- list()
  for (parm in names(TgaussParms$min)){
    set.seed(42)
    input_prosail[[parm]] <- truncnorm::rtruncnorm(n = nb_samples,
                                                   a = TgaussParms$min[[parm]],
                                                   b = TgaussParms$max[[parm]],
                                                   mean = TgaussParms$mean[[parm]],
                                                   sd = TgaussParms$sd[[parm]])
  }
  input_prosail <- data.frame(input_prosail)

  # define co-distribution with LAI
  if (codistribution_lai==TRUE){
    codistribution_lai <- list()
    codistribution_lai$Vmin0 <- data.frame('LIDFa' = 30, 'q' = 0.1, 'N' = 1.2,
                                           'CHL' = 20, 'LMA' = 0.003, 'Cw_rel' = 0.6,
                                           'BROWN' = 0.0, 'psoil' = 0)
    codistribution_lai$Vmax0 <- data.frame('LIDFa' = 80, 'q' = 0.5, 'N' = 1.8,
                                           'CHL' = 90, 'LMA' = 0.011, 'Cw_rel' = 0.85,
                                           'BROWN' = 2.0, 'psoil' = 1)
    codistribution_lai$VminLAImax <- data.frame('LIDFa' = 55, 'q' = 0.1, 'N' = 1.3,
                                                'CHL' = 45, 'LMA' = 0.005,
                                                'Cw_rel' = 0.70,'BROWN' = 0.0,
                                                'psoil' = 0)
    codistribution_lai$VmaxLAImax <- data.frame('LIDFa' = 65, 'q' = 0.5, 'N' = 1.8,
                                                'CHL' = 90, 'LMA' = 0.011,
                                                'Cw_rel' = 0.80, 'BROWN' = 0.2,
                                                'psoil' = 0.4)
    for (parm in names(codistribution_lai$Vmin0)){
      Vstar <- get_codistributions(V = input_prosail[[parm]],
                                   LAI = input_prosail$lai,
                                   MaxLAI = TgaussParms$max$lai,
                                   Vmin0 = codistribution_lai$Vmin0[[parm]],
                                   Vmax0 = codistribution_lai$Vmax0[[parm]],
                                   VminLAImax = codistribution_lai$VminLAImax[[parm]],
                                   VmaxLAImax = codistribution_lai$VmaxLAImax[[parm]])
      input_prosail[[parm]] <- Vstar
    }
  }

  # convert CW_rel into EWT
  input_prosail$EWT <- ((input_prosail$LMA)/(1-input_prosail$Cw_rel)) -
    input_prosail$LMA
  # set ANT, PROT and CBC to 0
  input_prosail$ANT <- input_prosail$PROT <- input_prosail$CBC <- 0
  # set CAR to 0.25*CHL
  input_prosail$CAR <- 0.25*input_prosail$CHL
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
  input_prosail$tts <- runif(n = nb_samples, min = GeomAcq['tts', 'min'],
                             max = GeomAcq['tts', 'max'])
  input_prosail$tto <- runif(n = nb_samples, min = GeomAcq['tto', 'min'],
                             max = GeomAcq['tto', 'max'])
  input_prosail$psi <- runif(n = nb_samples, min = GeomAcq['psi', 'min'],
                             max = GeomAcq['psi', 'max'])
  # default values
  input_prosail$TypeLidf <- 2
  input_prosail$alpha <- 40
  return(input_prosail)
}
