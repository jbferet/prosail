#' This function performs full training for hybrid inversion using SVR with
#' values for default parameters
#'
#' @param input_prosail list. user-defined list of input parameters to be used
#' to produce a training LUT
#' @param BRF_LUT list. user-defined BRF LUT used to run the hybrid inversion
#' @param atbd boolean. should input parameter distribution from ATBD be
#' applied ?
#' @param GeomAcq list. geometry of acquisition: min and max tts, tto & psi
#' @param codistribution_lai boolean. TRUE: codistrib with LAI accounted for
#' @param minval list. min val for input parms sampled to produce a training LUT
#' @param maxval list. max val for input parms sampled to produce a training LUT
#' @param TypeDistrib  list. Type of distribution. 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean value and STD corresponding to the
#' parameters sampled with gaussian distribution
#' @param parm_set list. list of input parameters set to a specific value
#' @param SAILversion character. Either 4SAIL or 4SAIL2
#' @param brown_lop character. Either 4SAIL or 4SAIL2
#' @param nb_samples numeric. number of samples in training LUT
#' @param nb_models numeric. number of individual models to be run for ensemble
#' @param replacement bolean. is there replacement in subsampling?
#' @param Parms2Estimate list. list of input parameters to be estimated
#' @param selected_bands list. bands used for regression for each input param
#' @param noise_level list. noise added to reflectance (defined per input parm)
#' @param SRF list. Spectral response function
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
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
#' input_variables based on BRF_LUT
#' @export

train_prosail_inversion <- function(input_prosail = NULL, BRF_LUT = NULL,
                                    atbd = FALSE, GeomAcq = NULL, SRF = NULL,
                                    codistribution_lai = TRUE, minval = NULL,
                                    maxval = NULL, TypeDistrib = NULL,
                                    GaussianDistrib = NULL, parm_set = NULL,
                                    SAILversion = '4SAIL', brown_lop = NULL,
                                    nb_samples = 2000, nb_models = 20,
                                    replacement = TRUE, Parms2Estimate = 'lai',
                                    selected_bands = NULL, noise_level = NULL,
                                    SpecPROSPECT = NULL, SpecSOIL = NULL,
                                    SpecATM = NULL, output_dir = './',
                                    method = 'liquidSVM', verbose = FALSE){

  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # define bands to select for inversion of each Parms2Estimate
  if (is.null(selected_bands)){
    Bands2Select <- NULL
  } else {
    Bands2Select <- list()
    if (is.null(SRF)){
      message('SRF not defined, cannot select bands according to their name')
      message('"selected_bands" and "SRF$Spectral_Bands" needs to partly match')
      message('prosail inversion will be trained with all variables available')
      Bands2Select <- NULL
    } else {
      for (parm in Parms2Estimate)
        Bands2Select[[parm]] <- match(selected_bands[[parm]], SRF$Spectral_Bands)
    }
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###     1- DEFINE THE LUT USED TO TRAIN THE HYBRID INVERSION          ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # default parameter values
  defaultVal <- data.frame('CHL'=40, 'CAR'=10, 'ANT' = 0, 'EWT' = 0.01,
                           'LMA' = 0.01, 'BROWN'=0.0, 'N' = 1.5, 'psoil' = 0.5,
                           'LIDFa' = 60, 'lai' = 2.5, 'q'=0.1, 'tto' = 0,
                           'tts' = 30, 'psi' = 80, 'TypeLidf' = 2, 'alpha' = 40)
  ListParms <- names(defaultVal)

  ##############################################################################
  #                     user-defined set of input parameters                  ##
  ##############################################################################
  if (!is.null(input_prosail)){
    input_prosail <- data.frame(input_prosail)
    if (atbd == TRUE | !is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('parms to generate BRF LUT provided by user in "input_prosail"')
        message('following input variables will be ignored: ')
        message('"atbd" "minval" "maxval" "TypeDistrib" "GaussianDistrib"')
      }
    }
    # check if all parameters defined & use default value for undefined parms
    UndefinedParms <- ListParms[which(is.na(match(ListParms,
                                                  names(input_prosail))))]
    for (parm in UndefinedParms) input_prosail[[parm]] <- defaultVal[[parm]]
    # Set parameters
    if (length(parm_set)>0){
      for (parm in names(parm_set)){
        if (is.null(input_prosail[[parm]]))
          input_prosail[[parm]] <- parm_set[[parm]]
      }
    }
  } else {
    input_prosail <- get_input_PROSAIL(atbd = atbd, GeomAcq = GeomAcq,
                                       codistribution_lai = codistribution_lai,
                                       minval = minval, maxval = maxval,
                                       TypeDistrib = TypeDistrib,
                                       GaussianDistrib = GaussianDistrib,
                                       parm_set = parm_set,
                                       nb_samples = nb_samples,
                                       verbose = verbose)
  }

  if (!is.null(BRF_LUT)){
    for (parm in Parms2Estimate){
      if (!is.null(BRF_LUT[[parm]])){
        BRF_LUT_Noise[[parm]] <- BRF_LUT[[parm]]
      } else {
        message('Please make sure you provide BRF_LUT as a list with elements corresponding to Parms2Estimate')
      }
    }
  } else {
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    ### 2- PRODUCE BRF from input_prosail & ddefault spectral sampling = 1nm  ##
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    # define default SpecPROSPECT, SpecSOIL and SpecATM if undefined
    if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
    if (is.null(SpecSOIL)) SpecSOIL <- prosail::SpecSOIL
    if (is.null(SpecATM)) SpecATM <- prosail::SpecATM
    # check if same spectral sampling for all key variables
    check_spectral_sampling(SpecPROSPECT, SpecSOIL, SpecATM)
    # generate LUT of BRF corresponding to input_prosail, for a sensor

    if (!'fCover' %in% Parms2Estimate &
        !'albedo' %in% Parms2Estimate &
        !'fAPAR' %in% Parms2Estimate){
      BRF_LUT <- generate_LUT_BRF(SAILversion = SAILversion,
                                  input_prosail = input_prosail,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM,
                                  brown_lop = brown_lop)
    } else if ('fCover' %in% Parms2Estimate |
               'albedo' %in% Parms2Estimate |
               'fAPAR' %in% Parms2Estimate){
      res <- generate_LUT_PROSAIL(SAILversion = SAILversion,
                                  input_prosail = input_prosail,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM,
                                  brown_lop = brown_lop)
      BRF_LUT <- res$BRF
      input_prosail$fCover <- res$fCover
      input_prosail$fAPAR <- res$fAPAR
      input_prosail$albedo <- res$albedo
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     3- APPLY SPECTRAL RESPONSE FUNCTION if not already applied    ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # apply sensor spectral response function if provided
    wvl <- SpecPROSPECT$lambda
    if (!is.null(SRF)) {
      if (!length(SRF$Spectral_Bands)==nrow(BRF_LUT)){
        BRF_LUT <- apply_sensor_characteristics(wvl = wvl, SRF = SRF,
                                                input_refl_table = BRF_LUT)
        SpecSensor <- prepare_sensor_simulation(SpecPROSPECT, SpecSOIL,
                                                SpecATM, SRF)
      }
      rownames(BRF_LUT) <- SRF$Spectral_Bands
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
    write.table(x = format(t(BRF_LUT), digits=5), file = filename,
                append = FALSE, quote = FALSE,
                col.names = SpecSensor$SpecPROSPECT_Sensor$lambda,
                row.names = FALSE,sep = '\t')

    # bands used for inversion
    for (parm in Parms2Estimate){
      if (is.null(Bands2Select[[parm]]))
        Bands2Select[[parm]] <- seq_len(length(SpecSensor$band_names))
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     4- add noise to reflectance data                              ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # if noise_level == NULL then use the same strategy than ATBD
    BRF_LUT_Noise <- list()
    if (is.null(noise_level)){
      if (SRF$Sensor %in% c('Sentinel_2', 'Sentinel_2A', 'Sentinel_2B')){
        BRF_LUT_NoiseAll <- apply_noise_atbd(BRF_LUT)
        for (parm in Parms2Estimate)
          BRF_LUT_Noise[[parm]] <- BRF_LUT_NoiseAll[Bands2Select[[parm]],]
      } else {
        for (parm in Parms2Estimate) {
          noise_level[[parm]] <- 0.01
          subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
          BRF_LUT_Noise[[parm]] <- subsetRefl +
            subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                    mean = 0,
                                    sd = noise_level[[parm]]),
                              nrow = nrow(subsetRefl))
        }
      }
    } else {
      # produce LUT with noise
      for (parm in Parms2Estimate){
        if (is.null(noise_level[[parm]])) noise_level[[parm]] <- 0.01
        subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
        BRF_LUT_Noise[[parm]] <- subsetRefl +
          subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                  mean = 0,
                                  sd = noise_level[[parm]]),
                            nrow = nrow(subsetRefl))
      }
    }
  }

  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     PERFORM HYBRID INVERSION                      ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # train SVR for each variable and each run
  modelSVR <- list()
  for (parm in Parms2Estimate){
    ColParm <- which(parm==names(input_prosail))
    input_variables <- input_prosail[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                             input_variables = input_variables,
                                             nb_bagg = nb_models,
                                             replacement = replacement,
                                             method = method, verbose = verbose)
  }
  return(modelSVR)
}
