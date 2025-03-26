test_that("PROSAIL iterative optimization", {
  # define lower and upper bounds for inversion, as well as initial value
  LB <- data.frame('CHL' = 5, 'CAR' = 1, 'EWT' = 0.002, 'LMA' = 0,
                   'lai' = 0.5, 'N' = 1)
  UB <- data.frame('CHL' = 80, 'CAR' = 20, 'EWT' = 0.03, 'LMA' = 0.03,
                   'lai' = 6, 'N' = 3)
  Init <- data.frame('CHL' = 40, 'CAR' = 10, 'EWT' = 0.01, 'LMA' = 0.01,
                     'lai' = 3, 'N' = 1.5)
  # define parameters which are already set for inversion
  ParmSet <- data.frame('tts' = 40, 'tto' = 0, 'psi' = 60,  'psoil' = 0,
                        'LIDFa' = 60, 'ANT' = 0, 'BROWN' = 0, 'q' = 0.1)
  # compute soil reflectance
  rsoil <- ParmSet$psoil*SpecSOIL$Dry_Soil+(1-ParmSet$psoil)*SpecSOIL$Wet_Soil
  # simulate canopy BRF with 1 nm sampling
  truth <- data.frame('CHL' = 60, 'CAR' = 8, 'EWT' = 0.015, 'LMA' = 0.005,
                      'lai' = 5,  'N' = 1.8)
  Refl_1nm <- PRO4SAIL(N = truth$N, CHL = truth$CHL, CAR = truth$CAR,
                       ANT = ParmSet$ANT, BROWN = ParmSet$BROWN, EWT = truth$EWT,
                       LMA = truth$LMA, TypeLidf = 2, lai = truth$lai,
                       q = ParmSet$q, LIDFa = ParmSet$LIDFa, rsoil = rsoil,
                       tts = ParmSet$tts, tto = ParmSet$tto, psi = ParmSet$psi)
  brf_1nm <- prosail::Compute_BRF(rdot = Refl_1nm$rdot,
                                  rsot = Refl_1nm$rsot,
                                  tts = ParmSet$tts,
                                  SpecATM_Sensor = SpecATM)
  # invert 1 nm data
  est <- Invert_PROSAIL(brfMES = brf_1nm$BRF,
                        InitialGuess = Init,
                        LowerBound = LB,
                        UpperBound = UB,
                        SpecPROSPECT_Sensor = SpecPROSPECT_FullRange,
                        SpecATM_Sensor = SpecATM,
                        SpecSOIL_Sensor = SpecSOIL,
                        TypeLidf = 2, ParmSet = ParmSet)

  nerr <- list()
  for (parm in names(truth))
    nerr[[parm]] <- 100*abs((truth[[parm]]-est[[parm]])/truth[[parm]])
  nerr <- unlist(nerr)
  expect_true(max(nerr)<1)
})
