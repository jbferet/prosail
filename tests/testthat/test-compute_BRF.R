test_that("prosail produces physically possible BRF values", {
  # run PROSAIL with 4SAIL2
  refl <- prosail()
  brf_4sail <- compute_brf(rdot = refl$rdot,
                           rsot = refl$rsot,
                           tts = 30,
                           spec_atm_sensor = spec_atm)

  # run PROSAIL with 4SAIL2
  input_prospect <- data.frame('chl' = c(40, 5), 'car' = c(8, 4),
                               'ant' = c(0.0, 1), 'ewt' = c(0.02, 0.01),
                               'lma' = c(0.009, 0.009), 'n_struct' = c(1.5, 2),
                               'brown' = c(0, 1))
  refl_4sail2 <- prosail(SAILversion = '4SAIL2', input_prospect= input_prospect,
                         type_lidf = 2, lidf_a = 30, lai = 5, hotspot = 0.1, tts = 30,
                         tto = 10, psi = 90, rsoil = spec_soil$max_refl,
                         fraction_brown = 0.5, diss = 0.5, cv = 1, zeta = 1)
  brf_4sail2 <- compute_brf(rdot = refl_4sail2$rdot,
                            rsot = refl_4sail2$rsot,
                            tts = 30,
                            spec_atm_sensor = spec_atm)

  expect_true(all(brf_4sail$BRF >= 0))
  expect_true(all(brf_4sail2$BRF >= 0))
})
