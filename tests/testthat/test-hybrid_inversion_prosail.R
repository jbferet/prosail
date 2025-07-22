test_that("hybrid inversion ok", {
  # get sensor response for Sentinel-2
  sensor_name <- 'Sentinel_2'
  srf <- get_radiometry(sensor_name)
  # define parameters to estimate
  Parms2Estimate <- 'lai'

  # define spectral bands required to train SVR model for each variable
  S2BandSelect <- c('B3','B4','B8')
  Bands2Select <- list('lai' = match(S2BandSelect,srf$spectral_bands))

  # define output directory where LUTs will be saved
  PROSAIL_ResPath <- './'

  # define ranges for geometry of acquisition
  geom_acq <- list()
  geom_acq$min <- data.frame('tto' = 0, 'tts' = 20, 'psi' = 0)
  geom_acq$max <- data.frame('tto' = 10, 'tts' = 30, 'psi' = 360)

  # define inputs for PROSAIL following ATBD
  input_prosail <- get_input_prosail(atbd = TRUE, geom_acq = geom_acq)

  # produce LUT
  res <- generate_lut_prosail(SAILversion = '4SAIL',
                              input_prosail = input_prosail,
                              spec_prospect = spec_prospect_fullrange,
                              spec_soil = prosail::spec_soil,
                              spec_atm = prosail::spec_atm)
  brf_lut_1nm <- res$brf
  brf_lut <- apply_sensor_characteristics(wvl = spec_prospect_fullrange$lambda,
                                          srf = srf,
                                          input_refl_table = brf_lut_1nm)
  # identify spectral bands in LUT
  rownames(brf_lut) <- srf$spectral_bands
  # add noise
  subset_refl <- brf_lut[Bands2Select$lai,]
  brf_lut_noise <- apply_noise_atbd(subset_refl)
  # train model
  modelSVR <- prosail_hybrid_train(brf_lut = brf_lut_noise,
                                   input_variables = input_prosail$lai,
                                   method = 'svmLinear')

  ###########################################################
  # perform prediction based on models in previous steps
  # the prediction returns mean value obtained form the ensemble of regression
  # models for each sample, as well as corresponding standard deviation
  hybrid_res <- prosail_hybrid_apply(regression_models = modelSVR,
                                    refl = brf_lut_noise)
  expect_true(cor.test(hybrid_res$MeanEstimate, input_prosail$lai)$estimate>0.7)
})
