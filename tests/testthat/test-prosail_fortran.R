test_that("compare prosail outputs with fortran distribution", {

  input_prospect <- list('chl' =	40, 'car' =	8, 'ant' =	0.5, 'brown' = 0.0,
                         'ewt' =	0.01, 'lma'  = 0.009, 'nstruct' =	1.5)
  lop <- prospect::prospect(input_prospect = input_prospect)
  rsoil <- spec_soil$max_refl
  rcan <- fourSAIL(lop = lop, type_lidf = 2, lidf_a = 30, lidf_b = NULL,
                   lai = 3, hotspot = 0.01, tts = 30, tto = 10, psi = 0,
                   rsoil = rsoil)
  rs <- get_surf_refl(rcan$rdot, rcan$rsot, tts = 30,
                      spec_atm_sensor = spec_atm, skyl = NULL)

  rf <- read.delim(file = './reflectance_fortran.txt' , sep = '\t')
  expect_true(sum(abs(rf$surf_refl- rs$surf_refl))<0.005)
})
