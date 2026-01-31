test_that("PROSAIL iterative optimization", {
  # define lower and upper bounds for inversion, as well as initial value
  lb <- data.frame('chl' = 5, 'car' = 1, 'ewt' = 0.002, 'lma' = 0,
                   'lai' = 0.5, 'n_struct' = 1)
  ub <- data.frame('chl' = 80, 'car' = 20, 'ewt' = 0.03, 'lma' = 0.03,
                   'lai' = 6, 'n_struct' = 3)
  init <- data.frame('chl' = 40, 'car' = 10, 'ewt' = 0.01, 'lma' = 0.01,
                     'lai' = 3, 'n_struct' = 1.5)
  # define parameters which are already set for inversion
  parm_set <- data.frame('tts' = 40, 'tto' = 0, 'psi' = 60,
                         'soil_brightness' = 1, 'lidf_a' = 60, 'ant' = 0,
                         'brown' = 0, 'hotspot' = 0.1)
  # compute soil reflectance
  rsoil <- parm_set$soil_brightness*spec_soil_ossl$soil_01
  # simulate canopy surf_refl with 1 nm sampling
  truth <- data.frame('chl' = 60, 'car' = 8, 'ewt' = 0.015, 'lma' = 0.005,
                      'lai' = 5,  'n_struct' = 1.8)
  refl_1nm <- prosail(n_struct = truth$n_struct, chl = truth$chl,
                      car = truth$car, ant = parm_set$ant,
                      brown = parm_set$brown, ewt = truth$ewt, lma = truth$lma,
                      type_lidf = 2, lai = truth$lai, hotspot = parm_set$hotspot,
                      lidf_a = parm_set$lidf_a, rsoil = rsoil,
                      tts = parm_set$tts, tto = parm_set$tto,
                      psi = parm_set$psi)
  surf_refl_1nm <- compute_surf_refl(rdot = refl_1nm$rdot,
                                     rsot = refl_1nm$rsot,
                                     tts = parm_set$tts,
                                     spec_atm_sensor = spec_atm)
  # invert 1 nm data
  spec_soil_sensor <- spec_soil_ossl[c('lambda', 'soil_01')]
  names(spec_soil_sensor) <- c('lambda', 'refl')
  est <- invert_prosail(refl_mes = surf_refl_1nm$surf_refl,
                        initialization = init,
                        lower_bound = lb,
                        upper_bound = ub,
                        spec_prospect_sensor = spec_prospect_full_range,
                        spec_atm_sensor = spec_atm,
                        spec_soil_sensor = spec_soil_sensor,
                        type_lidf = 2, parm_set = parm_set)

  nerr <- list()
  for (parm in names(truth))
    nerr[[parm]] <- 100*abs((truth[[parm]]-est[[parm]])/truth[[parm]])
  nerr <- unlist(nerr)
  expect_true(max(nerr)<1)
})
