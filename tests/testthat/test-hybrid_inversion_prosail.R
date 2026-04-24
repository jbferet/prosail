test_that("hybrid inversion ok", {
  # get sensor response for Sentinel-2
  srf <- get_srf_sensor(sensor_name = 'sentinel_2')
  # define spectral bands required to train SVR model for each variable
  s2_bands_to_select <- c('B3','B4','B8')
  bands_to_select <- list('lai' = match(s2_bands_to_select,srf$spectral_bands))
  # define output directory where LUTs will be saved
  PROSAIL_ResPath <- './'
  # define ranges for geometry of acquisition
  geom_acq <- list()
  geom_acq$min <- data.frame('tto' = 0, 'tts' = 20, 'psi' = 0)
  geom_acq$max <- data.frame('tto' = 10, 'tts' = 30, 'psi' = 360)
  # define inputs for PROSAIL following ATBD
  input_prosail <- get_input_prosail(atbd = TRUE, geom_acq = geom_acq)
  # produce LUT
  spec_prospect_full_range <- prosail::spec_prospect_full_range
  res <- generate_lut_prosail(SAILversion = '4SAIL',
                              input_prosail = input_prosail,
                              spec_prospect = spec_prospect_full_range,
                              spec_soil = prosail::spec_soil_atbd_v2,
                              spec_atm = prosail::spec_atm)
  refl_lut_1nm <- res$surf_refl
  refl_lut <- apply_sensor_characteristics(wvl = spec_prospect_full_range$lambda,
                                           srf = srf,
                                           refl = refl_lut_1nm)
  # identify spectral bands in LUT
  rownames(refl_lut) <- srf$spectral_bands
  # add noise
  subset_refl <- refl_lut[bands_to_select$lai,]
  refl_lut_noise <- apply_noise_atbd(refl_lut = subset_refl)
  # train model
  hybrid_model <- prosail_hybrid_train(refl_lut = refl_lut_noise,
                                       input_variables = input_prosail$lai,
                                       method = 'nu-svr')

  ###########################################################
  # perform prediction based on models in previous steps
  # the prediction returns mean value obtained form the ensemble of regression
  # models for each sample, as well as corresponding standard deviation
  hybrid_res <- prosail_hybrid_apply(regression_models = hybrid_model,
                                     refl = refl_lut_noise)
  expect_true(cor.test(hybrid_res$mean_estimate, input_prosail$lai)$estimate>0.7)
})
