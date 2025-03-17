test_that("PRO4AIL and PRO4SAIL2 produce physically possible BRF values", {
  # run PROSAIL with 4SAIL2
  Refl <- PRO4SAIL()
  BRF_4SAIL <- compute_BRF(rdot = Refl$rdot,
                           rsot = Refl$rsot,
                           tts = 30,
                           SpecATM_Sensor = prosail::SpecATM)

  # run PROSAIL with 4SAIL2
  Input_PROSPECT <- data.frame('CHL' = c(40, 5), 'CAR' = c(8, 4),
                               'ANT' = c(0.0, 1), 'EWT' = c(0.02, 0.01),
                               'LMA' = c(0.009, 0.009), 'N' = c(1.5, 2),
                               'BROWN' = c(0, 1))
  Ref_4SAIL2 <- PRO4SAIL(SAILversion = '4SAIL2', Input_PROSPECT= Input_PROSPECT,
                         TypeLidf = 2, LIDFa = 30, lai = 5, q = 0.1, tts = 30,
                         tto = 10, psi = 90, rsoil = prosail::SpecSOIL$Dry_Soil,
                         fraction_brown = 0.5, diss = 0.5, Cv = 1, Zeta = 1)
  BRF_4SAIL2 <- compute_BRF(rdot = Ref_4SAIL2$rdot,
                            rsot = Ref_4SAIL2$rsot,
                            tts = 30,
                            SpecATM_Sensor = prosail::SpecATM)

  expect_true(all(BRF_4SAIL$BRF >= 0))
  expect_true(all(BRF_4SAIL2$BRF >= 0))
})
