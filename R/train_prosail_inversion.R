#' This function performs full training for hybrid inversion using SVR with
#' values for default parameters
#'
#' @param InputPROSAIL list. user-defined list of input parameters to be used
#' to produce a training LUT
#' @param BRF_LUT list. user-defined BRF LUT used to run the hybrid inversion
#' @param atbd boolean. should input parameter distribution from ATBD be
#' applied ?
#' @param GeomAcq list. geometry of acquisition: min and max tts, tto & psi
#' @param Codist_LAI boolean. set TYRUE if codistribution with LAI accounted for
#' @param minval list. min val for input parms sampled to produce a training LUT
#' @param maxval list. max val for input parms sampled to produce a training LUT
#' @param TypeDistrib  list. Type of distribution. 'Uniform' or 'Gaussian'
#' @param GaussianDistrib  list. Mean value and STD corresponding to the
#' parameters sampled with gaussian distribution
#' @param ParmSet list. list of input parameters set to a specific value
#' @param SAILversion character. Either 4SAIL or 4SAIL2
#' @param BrownLOP character. Either 4SAIL or 4SAIL2
#' @param nbSamples numeric. number of samples in training LUT
#' @param nbSamplesPerRun numeric. nb of training sample per individual regression model
#' @param nbModels numeric. number of individual models to be run for ensemble
#' @param Replacement bolean. is there replacement in subsampling?
#' @param Parms2Estimate list. list of input parameters to be estimated
#' @param Bands2Select list. bands used for regression for each input parameter
#' @param NoiseLevel list. noise added to reflectance (defined per input parm)
#' @param SRF list. Spectral response function
#' @param SpecPROSPECT list. Includes optical constants required for PROSPECT
#' @param SpecSOIL list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param SpecATM list. direct and diffuse radiation for clear conditions
#' @param Path_Results character. path for results
#' @param FigPlot boolean. Set TRUE to get scatterplot of estimated biophysical
#' variable during training step
#' @param method character. which machine learning regression method should be
#' used?
#' default = SVM with liquidSVM. svmRadial and svmLinear from caret package also
#' implemented. More to come
#' @param verbose boolean. when set to TRUE, prints message if hyperparameter
#' adjustment performed during training
#'
#' @return modelsSVR list. regression models trained for the retrieval of
#' InputVar based on BRF_LUT
#' @export

train_prosail_inversion <- function(InputPROSAIL = NULL, BRF_LUT = NULL,
                                    atbd = FALSE, GeomAcq = NULL,
                                    Codist_LAI = TRUE, minval = NULL,
                                    maxval = NULL, TypeDistrib = NULL,
                                    GaussianDistrib = NULL, ParmSet = NULL,
                                    SAILversion = '4SAIL', BrownLOP = NULL,
                                    nbSamples = 2000, nbSamplesPerRun = 100,
                                    nbModels = 20, Replacement = TRUE,
                                    Parms2Estimate = 'lai', Bands2Select = NULL,
                                    NoiseLevel = NULL, SRF = NULL,
                                    SpecPROSPECT = NULL, SpecSOIL = NULL,
                                    SpecATM = NULL, Path_Results = './',
                                    FigPlot = FALSE, method = 'liquidSVM',
                                    verbose = FALSE){

  dir.create(Path_Results, showWarnings = FALSE, recursive = TRUE)
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
  if (!is.null(InputPROSAIL)){
    InputPROSAIL <- data.frame(InputPROSAIL)
    if (atbd == TRUE | !is.null(minval) | !is.null(maxval)){
      if (verbose==TRUE){
        message('parameters to generate BRF LUT provided by user in "InputPROSAIL"')
        message('following input variables will be ignored: "atbd" "minval" "maxval" "TypeDistrib" "GaussianDistrib"')
      }
    }
    # check if all parameters defined & use default value for undefined parameters
    UndefinedParms <- ListParms[which(is.na(match(ListParms,
                                                  names(InputPROSAIL))))]
    for (parm in UndefinedParms) InputPROSAIL[[parm]] <- defaultVal[[parm]]
    # Set parameters
    if (length(ParmSet)>0){
      for (parm in names(ParmSet)){
        if (is.null(InputPROSAIL[[parm]]))
          InputPROSAIL[[parm]] <- ParmSet[[parm]]
      }
    }
  } else {
    InputPROSAIL <- get_input_PROSAIL(atbd = atbd, GeomAcq = GeomAcq,
                                     Codist_LAI = Codist_LAI,
                                     minval = minval, maxval = maxval,
                                     TypeDistrib = TypeDistrib,
                                     GaussianDistrib = GaussianDistrib,
                                     ParmSet = ParmSet, nbSamples = nbSamples,
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
    ### 2- PRODUCE BRF from InputPROSAIL & ddefault spectral sampling = 1nm  ###
    ### == == == == == == == == == == == == == == == == == == == == == == == ###
    # define default SpecPROSPECT, SpecSOIL and SpecATM if undefined
    if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
    if (is.null(SpecSOIL)) SpecSOIL <- prosail::SpecSOIL
    if (is.null(SpecATM)) SpecATM <- prosail::SpecATM
    # check if same spectral sampling for all key variables
    check_spectral_sampling(SpecPROSPECT, SpecSOIL, SpecATM)
    # generate LUT of BRF corresponding to InputPROSAIL, for a sensor

    if (!'fCover' %in% Parms2Estimate &
        !'albedo' %in% Parms2Estimate &
        !'fAPAR' %in% Parms2Estimate){
      BRF_LUT <- generate_LUT_BRF(SAILversion = SAILversion,
                                  InputPROSAIL = InputPROSAIL,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM,
                                  BrownLOP = BrownLOP)
    } else if ('fCover' %in% Parms2Estimate |
               'albedo' %in% Parms2Estimate |
               'fAPAR' %in% Parms2Estimate){
      res <- generate_LUT_PROSAIL(SAILversion = SAILversion,
                                  InputPROSAIL = InputPROSAIL,
                                  SpecPROSPECT = SpecPROSPECT,
                                  SpecSOIL = SpecSOIL,
                                  SpecATM = SpecATM,
                                  BrownLOP = BrownLOP)
      BRF_LUT <- res$BRF
      InputPROSAIL$fCover <- res$fCover
      InputPROSAIL$fAPAR <- res$fAPAR
      InputPROSAIL$albedo <- res$albedo
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     3- APPLY SPECTRAL RESPONSE FUNCTION if not already applied    ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # apply sensor spectral response function if provided
    wvl <- SpecPROSPECT$lambda
    if (!is.null(SRF)) {
      if (!length(SRF$Spectral_Bands)==nrow(BRF_LUT)){
        BRF_LUT <- apply_sensor_characteristics(wvl = wvl, SRF = SRF,
                                              InRefl = BRF_LUT)
        SpecSensor <- prepare_sensor_simulation(SpecPROSPECT, SpecSOIL,
                                              SpecATM, SRF)
      }
      rownames(BRF_LUT) <- SRF$Spectral_Bands
    }

    # write parameters LUT
    output <- matrix(unlist(InputPROSAIL), ncol = length(InputPROSAIL),
                     byrow = FALSE)
    filename <- file.path(Path_Results,'PROSAIL_LUT_InputParms.txt')
    write.table(x = format(output, digits=3), file = filename,append = FALSE,
                quote = FALSE, col.names = names(InputPROSAIL),
                row.names = FALSE, sep = '\t')
    # Write BRF LUT corresponding to parameters LUT
    filename <- file.path(Path_Results,'PROSAIL_LUT_Reflectance.txt')
    write.table(x = format(t(BRF_LUT), digits=5), file = filename,
                append = FALSE, quote = FALSE,
                col.names = SpecSensor$SpecPROSPECT_Sensor$lambda,
                row.names = FALSE,sep = '\t')

    # bands used for inversion
    for (parm in Parms2Estimate){
      if (is.null(Bands2Select[[parm]]))
        Bands2Select[[parm]] <- seq_len(length(SpecSensor$BandNames))
    }

    ### == == == == == == == == == == == == == == == == == == == == == == ###
    ###     4- add noise to reflectance data                              ###
    ### == == == == == == == == == == == == == == == == == == == == == == ###
    # if NoiseLevel == NULL then use the same strategy than ATBD
    BRF_LUT_Noise <- list()
    if (is.null(NoiseLevel)){
      if (SRF$Sensor %in% c('Sentinel_2', 'Sentinel_2A', 'Sentinel_2B')){
        BRF_LUT_NoiseAll <- apply_noise_atbd(BRF_LUT)
        for (parm in Parms2Estimate)
          BRF_LUT_Noise[[parm]] <- BRF_LUT_NoiseAll[Bands2Select[[parm]],]
      } else {
        for (parm in Parms2Estimate) {
          NoiseLevel[[parm]] <- 0.01
          subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
          BRF_LUT_Noise[[parm]] <- subsetRefl +
            subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                    mean = 0,
                                    sd = NoiseLevel[[parm]]),
                              nrow = nrow(subsetRefl))
        }
      }
    } else {
      # produce LUT with noise
      for (parm in Parms2Estimate){
        if (is.null(NoiseLevel[[parm]])) NoiseLevel[[parm]] <- 0.01
        subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
        BRF_LUT_Noise[[parm]] <- subsetRefl +
          subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),
                                  mean = 0,
                                  sd = NoiseLevel[[parm]]),
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
    ColParm <- which(parm==names(InputPROSAIL))
    InputVar <- InputPROSAIL[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                             InputVar = InputVar,
                                             nbEnsemble = nbModels,
                                             WithReplacement = Replacement,
                                             method = method, verbose = verbose)
  }
  return(modelSVR)
}
