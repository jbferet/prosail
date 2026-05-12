#' This function provides a LUT corresponding to PROSAIL input parameters to
#' reproduce the distribution of parameters used in the ATBD
#'
#' @param nb_samples numeric. number of samples to generate
#' @param geom_acq list. min and max values for tts, tto and psi
#' @param codistribution_lai boolean. set TYRUE if codistribution with LAI accounted for
#'
#' @return input_prosail list. list of PROSAIL input parameters
#' @importFrom truncnorm rtruncnorm
#' @export

get_atbd_lut_input <- function(nb_samples = 2000, geom_acq = NULL,
                               codistribution_lai = TRUE){
  # define paremertization for truncated gaussians
  t_gauss_parms <- list()
  t_gauss_parms$min <- data.frame('lai' = 0, 'lidf_a' = 30, 'hotspot' = 0.1,
                                  'n_struct' = 1.2, 'chl' = 20, 'lma' = 0.003,
                                  'cw_rel' = 0.6, 'brown' = 0.0, 'psoil' = 0,
                                  'soil_brightness' = 0.5)
  t_gauss_parms$max <- data.frame('lai' = 15, 'lidf_a' = 80, 'hotspot' = 0.5,
                                  'n_struct' = 1.8, 'chl' = 90, 'lma' = 0.011,
                                  'cw_rel' = 0.85, 'brown' = 2.0, 'psoil' = 1,
                                  'soil_brightness' = 3.5)
  t_gauss_parms$mean <- data.frame('lai' = 2, 'lidf_a' = 60, 'hotspot' = 0.2,
                                   'n_struct' = 1.5, 'chl' = 45, 'lma' = 0.005,
                                   # 'cw_rel' = 0.75, 'brown' = 0, 'psoil' = 0.25)
                                   'cw_rel' = 0.75, 'brown' = 0, 'psoil' = 0.5,
                                   'soil_brightness' = 1.2)
  t_gauss_parms$sd <- data.frame('lai' = 3, 'lidf_a' = 30, 'hotspot' = 0.5,
                                 'n_struct' = 0.3, 'chl' = 30, 'lma' = 0.005,
                                 'cw_rel' = 0.08, 'brown' = 0.30, 'psoil' = 0.6,
                                 'soil_brightness' = 2)

  # get distribution corresponding to gaussians
  input_prosail <- list()
  for (parm in names(t_gauss_parms$min)){
    set.seed(42)
    input_prosail[[parm]] <- truncnorm::rtruncnorm(n = nb_samples,
                                                   a = t_gauss_parms$min[[parm]],
                                                   b = t_gauss_parms$max[[parm]],
                                                   mean = t_gauss_parms$mean[[parm]],
                                                   sd = t_gauss_parms$sd[[parm]])
  }
  input_prosail <- data.frame(input_prosail)

  # define co-distribution with LAI
  if (codistribution_lai==TRUE){
    codistribution_lai <- list()
    codistribution_lai$vmin_0 <- data.frame('lidf_a' = 30, 'hotspot' = 0.1,
                                            'n_struct' = 1.2, 'chl' = 20,
                                            'lma' = 0.003, 'cw_rel' = 0.6,
                                            'brown' = 0.0, 'psoil' = 0,
                                            'soil_brightness' = 0.5)
    codistribution_lai$vmax_0 <- data.frame('lidf_a' = 80, 'hotspot' = 0.5,
                                            'n_struct' = 1.8, 'chl' = 90,
                                            'lma' = 0.011, 'cw_rel' = 0.85,
                                            'brown' = 2.0, 'psoil' = 1,
                                            'soil_brightness' = 3.5)
    codistribution_lai$vmin_lai_max <- data.frame('lidf_a' = 55, 'hotspot' = 0.1,
                                                  'n_struct' = 1.3, 'chl' = 45,
                                                  'lma' = 0.005, 'cw_rel' = 0.70,
                                                  'brown' = 0.0, 'psoil' = 0,
                                                  'soil_brightness' = 0.5)
    codistribution_lai$vmax_lai_max <- data.frame('lidf_a' = 65, 'hotspot' = 0.5,
                                                  'n_struct' = 1.8, 'chl' = 90,
                                                  'lma' = 0.011, 'cw_rel' = 0.80,
                                                  'brown' = 0.2, 'psoil' = 0.4,
                                                  'soil_brightness' = 1.2)
    for (parm in names(codistribution_lai$vmin_0)){
      v_star <- get_codistributions(V = input_prosail[[parm]],
                                    lai = input_prosail$lai,
                                    max_lai = t_gauss_parms$max$lai,
                                    vmin_0 = codistribution_lai$vmin_0[[parm]],
                                    vmax_0 = codistribution_lai$vmax_0[[parm]],
                                    vmin_lai_max = codistribution_lai$vmin_lai_max[[parm]],
                                    vmax_lai_max = codistribution_lai$vmax_lai_max[[parm]])
      input_prosail[[parm]] <- v_star
    }
  }

  # convert cw_rel into ewt
  input_prosail$ewt <- ((input_prosail$lma)/(1-input_prosail$cw_rel)) -
    input_prosail$lma
  # set ANT, PROT and CBC to 0
  input_prosail$ant <- input_prosail$prot <- input_prosail$cbc <- 0
  # set car to 0.25*chl
  input_prosail$car <- 0.25*input_prosail$chl
  # set geometry of acquisition
  if (is.null(geom_acq))
    geom_acq <- data.frame('min' = c('tto' = 0, 'tts' = 20, 'psi' = 0),
                           'max' = c('tto' = 10, 'tts' = 30, 'psi' = 360))
  if (inherits(geom_acq, "list"))
    geom_acq <- data.frame('min' = c('tto' = geom_acq$min$tto,
                                     'tts' = geom_acq$min$tts,
                                     'psi' = geom_acq$min$psi),
                           'max' = c('tto' = geom_acq$max$tto,
                                     'tts' = geom_acq$max$tts,
                                     'psi' = geom_acq$max$psi))
  input_prosail$tts <- runif(n = nb_samples, min = geom_acq['tts', 'min'],
                             max = geom_acq['tts', 'max'])
  input_prosail$tto <- runif(n = nb_samples, min = geom_acq['tto', 'min'],
                             max = geom_acq['tto', 'max'])
  input_prosail$psi <- runif(n = nb_samples, min = geom_acq['psi', 'min'],
                             max = geom_acq['psi', 'max'])

  # soil ID
  # input_prosail$soil_brightness <- 1
  input_prosail$soil_ID <- sample(x = 7, size = nb_samples, replace = T)

  # default values
  input_prosail$type_lidf <- 2
  input_prosail$alpha <- 40
  return(input_prosail)
}


#' @rdname prosail-deprecated
#' @export
get_atbd_LUT_input <- function(nbSamples = 2000, GeomAcq = NULL,
                               Codist_LAI = TRUE){
  .Deprecated("get_atbd_lut_input")
  get_atbd_lut_input(nbSamples, GeomAcq, Codist_LAI)
}
