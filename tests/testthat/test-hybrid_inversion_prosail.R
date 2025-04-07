test_that("hybrid inversion ok", {
  # get sensor response for Sentinel-2
  sensor_name <- 'Sentinel_2'
  SRF <- get_radiometry(sensor_name)
  # define parameters to estimate
  Parms2Estimate <- 'lai'

  # define spectral bands required to train SVR model for each variable
  S2BandSelect <- c('B3','B4','B8')
  Bands2Select <- list('lai' = match(S2BandSelect,SRF$Spectral_Bands))

  # define output directory where LUTs will be saved
  PROSAIL_ResPath <- './'

  # define ranges for geometry of acquisition
  GeomAcq <- list()
  GeomAcq$min <- data.frame('tto' = 0, 'tts' = 20, 'psi' = 0)
  GeomAcq$max <- data.frame('tto' = 10, 'tts' = 30, 'psi' = 360)

  # define inputs for PROSAIL following ATBD
  input_prosail <- get_input_PROSAIL(atbd = TRUE, GeomAcq = GeomAcq)

  # produce LUT
  res <- generate_LUT_PROSAIL(SAILversion = '4SAIL',
                              input_prosail = input_prosail,
                              SpecPROSPECT = SpecPROSPECT_FullRange,
                              SpecSOIL = prosail::SpecSOIL,
                              SpecATM = prosail::SpecATM)
  BRF_LUT_1nm <- res$BRF
  BRF_LUT <- apply_sensor_characteristics(wvl = SpecPROSPECT_FullRange$lambda,
                                          SRF = SRF,
                                          input_refl_table = BRF_LUT_1nm)
  # identify spectral bands in LUT
  rownames(BRF_LUT) <- SRF$Spectral_Bands
  # add noise
  subsetRefl <- BRF_LUT[Bands2Select$lai,]
  BRF_LUT_Noise <- apply_noise_atbd(subsetRefl)
  # train model
  modelSVR <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise,
                                   input_variables = input_prosail$lai,
                                   method = 'svmLinear')

  ###########################################################
  # perform prediction based on models in previous steps
  # the prediction returns mean value obtained form the ensemble of regression
  # models for each sample, as well as corresponding standard deviation
  HybridRes <- PROSAIL_Hybrid_Apply(RegressionModels = modelSVR,
                                    Refl = BRF_LUT_Noise)
  expect_true(cor.test(HybridRes$MeanEstimate, input_prosail$lai)$estimate>0.7)
})
