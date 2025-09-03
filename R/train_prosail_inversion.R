#' This function performs full training for hybrid inversion using SVR with
#' values for default parameters
#'
#' @param input_prosail list. user-defined list of input parameters to be used
#' to produce a training LUT
#' @param brf_lut list. user-defined BRF LUT used to run the hybrid inversion
#' @param atbd boolean. should input parameter distribution from ATBD be
#' applied ?
#' @param geom_acq list. geometry of acquisition: min and max tts, tto & psi
#' @param codistribution_lai boolean. TRUE: codistrib with LAI accounted for
#' @param minval list. min val for input parms sampled to produce a training LUT
#' @param maxval list. max val for input parms sampled to produce a training LUT
#' @param type_distrib  list. Type of distribution. 'Uniform' or 'Gaussian'
#' @param gaussian_distrib  list. Mean value and STD corresponding to the
#' parameters sampled with gaussian distribution
#' @param parm_set list. list of input parameters set to a specific value
#' @param SAILversion character. Either 4SAIL or 4SAIL2
#' @param brown_lop character. Either 4SAIL or 4SAIL2
#' @param nb_samples numeric. number of samples in training LUT
#' @param nb_models numeric. number of individual models to be run for ensemble
#' @param replacement bolean. is there replacement in subsampling?
#' @param parms_to_estimate list. list of input parameters to be estimated
#' @param selected_bands list. bands used for regression for each input param
#' @param noise_level list. noise added to reflectance (defined per input parm)
#' @param srf list. Spectral response function
#' @param spec_prospect list. Includes optical constants required for PROSPECT
#' @param spec_soil list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param spec_atm list. direct and diffuse radiation for clear conditions
#' @param output_dir character. path for results
#' variable during training step
#' @param method character. which machine learning regression method should be
#' used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also
#' implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#'
#' @return modelsSVR list. regression models trained for the retrieval of
#' input_variables based on brf_lut
#' @export

train_prosail_inversion <- function(input_prosail = NULL, brf_lut = NULL,
                                    atbd = FALSE, geom_acq = NULL, srf = NULL,
                                    codistribution_lai = TRUE, minval = NULL,
                                    maxval = NULL, type_distrib = NULL,
                                    gaussian_distrib = NULL, parm_set = NULL,
                                    SAILversion = '4SAIL', brown_lop = NULL,
                                    nb_samples = 2000, nb_models = 20,
                                    replacement = TRUE,
                                    parms_to_estimate = 'lai',
                                    selected_bands = NULL, noise_level = NULL,
                                    spec_prospect = NULL, spec_soil = NULL,
                                    spec_atm = NULL, output_dir = './',
                                    method = 'liquidSVM', verbose = FALSE){

  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # define bands to select for inversion of each parms_to_estimate
  if (is.null(selected_bands)){
    bands_to_select <- NULL
  } else {
    bands_to_select <- list()
    if (is.null(srf)){
      message('srf not defined, cannot select bands according to their name')
      message('"selected_bands" and "srf$spectral_bands" needs to partly match')
      message('prosail inversion will be trained with all variables available')
      bands_to_select <- NULL
    } else {
      for (parm in parms_to_estimate)
        bands_to_select[[parm]] <- match(selected_bands[[parm]],
                                         srf$spectral_bands)
    }
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###     1- DEFINE THE LUT USED TO TRAIN THE HYBRID INVERSION          ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # default parameter values
  default_val <- data.frame('chl'=40, 'car'=10, 'ant' = 0, 'ewt' = 0.01,
                            'lma' = 0.01, 'brown'=0.0, 'n_struct' = 1.5,
                            'psoil' = 0.5, 'lidf_a' = 60, 'lai' = 2.5, 'hotspot'=0.1,
                            'tto' = 0, 'tts' = 30, 'psi' = 80, 'type_lidf' = 2,
                            'alpha' = 40)
  list_parms <- names(default_val)

  ##############################################################################
  #                     user-defined set of input parameters                  ##
  ##############################################################################
  if (!is.null(input_prosail)){
    input_prosail <- data.frame(input_prosail)
    # ensure compatibility with variable names from previous versions
    input_prosail <- set_compatibility_input_prosail(input_prosail)
    if (atbd == TRUE | !is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('parms to generate BRF LUT provided by user in "input_prosail"')
        message('following input variables will be ignored: ')
        message('"atbd" "minval" "maxval" "type_distrib" "gaussian_distrib"')
      }
    }
    # check if all parameters defined & use default value for undefined parms
    undefined_parms <- list_parms[which(is.na(match(list_parms,
                                                    names(input_prosail))))]
    for (parm in undefined_parms) input_prosail[[parm]] <- default_val[[parm]]
    # Set parameters
    if (length(parm_set)>0){
      for (parm in names(parm_set)){
        if (is.null(input_prosail[[parm]]))
          input_prosail[[parm]] <- parm_set[[parm]]
      }
    }
  } else {
    input_prosail <- get_input_prosail(atbd = atbd, geom_acq = geom_acq,
                                       codistribution_lai = codistribution_lai,
                                       minval = minval, maxval = maxval,
                                       type_distrib = type_distrib,
                                       gaussian_distrib = gaussian_distrib,
                                       parm_set = parm_set,
                                       nb_samples = nb_samples,
                                       verbose = verbose)
  }

  if (!is.null(brf_lut)){
    for (parm in parms_to_estimate){
      if (!is.null(brf_lut[[parm]])){
        brf_lut_noise[[parm]] <- brf_lut[[parm]]
      } else {
        message('Please make sure brf_lut is a list with elements ')
        message('corresponding to parms_to_estimate')
      }
    }
  } else {
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    ### 2- PRODUCE BRF from input_prosail & ddefault spectral sampling = 1nm  ##
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    # define default spec_prospect, spec_soil and spec_atm if undefined
    if (is.null(spec_prospect))
      spec_prospect <- prosail::spec_prospect_fullrange
    if (is.null(spec_soil))
      spec_soil <- prosail::spec_soil
    if (is.null(spec_atm))
      spec_atm <- prosail::spec_atm
    # check if same spectral sampling for all key variables
    check_spectral_sampling(spec_prospect, spec_soil, spec_atm)
    # generate LUT of BRF corresponding to input_prosail, for a sensor

    if (!'fcover' %in% parms_to_estimate &
        !'albedo' %in% parms_to_estimate &
        !'fapar' %in% parms_to_estimate){
      brf_lut <- generate_lut_brf(SAILversion = SAILversion,
                                  input_prosail = input_prosail,
                                  spec_prospect = spec_prospect,
                                  spec_soil = spec_soil,
                                  spec_atm = spec_atm,
                                  brown_lop = brown_lop)
    } else if ('fcover' %in% parms_to_estimate |
               'albedo' %in% parms_to_estimate |
               'fapar' %in% parms_to_estimate){
      res <- generate_lut_prosail(SAILversion = SAILversion,
                                  input_prosail = input_prosail,
                                  spec_prospect = spec_prospect,
                                  spec_soil = spec_soil,
                                  spec_atm = spec_atm,
                                  brown_lop = brown_lop)
      brf_lut <- res$brf
      input_prosail$fcover <- res$fcover
      input_prosail$fapar <- res$fapar
      input_prosail$albedo <- res$albedo
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     3- APPLY SPECTRAL RESPONSE FUNCTION if not already applied    ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # apply sensor spectral response function if provided
    wvl <- spec_prospect$lambda
    if (!is.null(srf)) {
      if (!length(srf$spectral_bands)==nrow(brf_lut)){
        brf_lut <- apply_sensor_characteristics(wvl = wvl, srf = srf,
                                                input_refl_table = brf_lut)
        spec_sensor <- prepare_sensor_simulation(spec_prospect, spec_soil,
                                                 spec_atm, srf)
      }
      rownames(brf_lut) <- srf$spectral_bands
    }

    # write parameters LUT
    output <- matrix(unlist(input_prosail), ncol = length(input_prosail),
                     byrow = FALSE)
    filename <- file.path(output_dir,'PROSAIL_LUT_InputParms.txt')
    write.table(x = format(output, digits=3), file = filename,append = FALSE,
                quote = FALSE, col.names = names(input_prosail),
                row.names = FALSE, sep = '\t')
    # Write BRF LUT corresponding to parameters LUT
    filename <- file.path(output_dir,'PROSAIL_LUT_Reflectance.txt')
    write.table(x = format(t(brf_lut), digits=5), file = filename,
                append = FALSE, quote = FALSE,
                col.names = spec_sensor$spec_prospect_sensor$lambda,
                row.names = FALSE,sep = '\t')

    # bands used for inversion
    for (parm in parms_to_estimate){
      if (is.null(bands_to_select[[parm]]))
        bands_to_select[[parm]] <- seq_len(length(spec_sensor$band_names))
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     4- add noise to reflectance data                              ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # if noise_level == NULL then use the same strategy than ATBD
    brf_lut_noise <- list()
    if (is.null(noise_level)){
      if (srf$sensor %in% c('Sentinel_2', 'Sentinel_2A', 'Sentinel_2B')){
        brf_lut_noise_all <- apply_noise_atbd(brf_lut)
        for (parm in parms_to_estimate)
          brf_lut_noise[[parm]] <- brf_lut_noise_all[bands_to_select[[parm]],]
      } else {
        for (parm in parms_to_estimate) {
          noise_level[[parm]] <- 0.01
          refl_subset <- brf_lut[bands_to_select[[parm]],]
          brf_lut_noise[[parm]] <- refl_subset +
            refl_subset*matrix(rnorm(nrow(refl_subset)*ncol(refl_subset),
                                     mean = 0,
                                     sd = noise_level[[parm]]),
                               nrow = nrow(refl_subset))
        }
      }
    } else {
      # produce LUT with noise
      for (parm in parms_to_estimate){
        if (is.null(noise_level[[parm]])) noise_level[[parm]] <- 0.01
        refl_subset <- brf_lut[bands_to_select[[parm]],]
        brf_lut_noise[[parm]] <- refl_subset +
          refl_subset*matrix(rnorm(nrow(refl_subset)*ncol(refl_subset),
                                   mean = 0,
                                   sd = noise_level[[parm]]),
                             nrow = nrow(refl_subset))
      }
    }
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     PERFORM HYBRID INVERSION                      ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # train SVR for each variable and each run
  model_svr <- list()
  for (parm in parms_to_estimate){
    ColParm <- which(parm==names(input_prosail))
    input_variables <- input_prosail[[ColParm]]
    model_svr[[parm]] <- prosail_hybrid_train(brf_lut = brf_lut_noise[[parm]],
                                              input_variables = input_variables,
                                              nb_bagg = nb_models,
                                              replacement = replacement,
                                              method = method,
                                              verbose = verbose)
  }
  return(model_svr)
}
