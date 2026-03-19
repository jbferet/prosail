#' adjust prospect inputs and run prospect in preparation for 4SAIL or 4SAIL2
#'
#' @param sail_version character. should be '4SAIL' or '4SAIL2'
#' @param spec_sensor dataframe. spectral properties required to run PROSPECT
#' @param input_prospect dataframe. includes all prospect input parameters
#' @param chl numeric. Chlorophyll content (microg.cm-2)
#' @param car numeric. Carotenoid content (microg.cm-2)
#' @param ant numeric. Anthocyain content (microg.cm-2)
#' @param brown numeric. Brown pigment content (Arbitrary units)
#' @param ewt numeric. Equivalent Water Thickness (g.cm-2)
#' @param lma numeric. Leaf Mass per Area (g.cm-2)
#' @param prot numeric. protein content  (g.cm-2)
#' @param cbc numeric. Carbon-based constituent content (g.cm-2)
#' @param n_struct numeric. Leaf structure parameter
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' (simulation of roughness)
#' @param fraction_brown numeric. fraction of brown vegetation (0-1)
#' @param brown_lop dataframe. brown leaf optical properties, when available
#'
#' @return invisible
#' @export

adjust_prospect_to_sail <- function(sail_version, spec_sensor, input_prospect,
                                    chl, car, ant, brown, ewt, lma, prot, cbc,
                                    n_struct, alpha, fraction_brown,
                                    brown_lop = NULL){

  # for all versions of 4SAIL: get green vegetation
  if (utils::packageVersion("prospect")<'2.0.0'){
    inprospect_green <- prospect::define_Input_PROSPECT(input_prospect[1,],
                                                        chl[1], car[1], ant[1],
                                                        brown[1], ewt[1], lma[1],
                                                        prot[1], cbc[1],
                                                        n_struct[1], alpha[1])
    green_lop <- prospect::PROSPECT(SpecPROSPECT = spec_sensor,
                                    Input_PROSPECT = inprospect_green)
    names(green_lop) <- c('wvl', 'reflectance', 'transmittance')
  } else {
    inprospect_green <- prospect::define_input_prospect(input_prospect = input_prospect[1,],
                                                        chl = chl[1],
                                                        car = car[1],
                                                        ant = ant[1],
                                                        brown = brown[1],
                                                        ewt = ewt[1],
                                                        lma = lma[1],
                                                        prot = prot[1],
                                                        cbc = cbc[1],
                                                        n_struct = n_struct[1],
                                                        alpha = alpha[1])
    green_lop <- prospect::prospect(spec_prospect = spec_sensor,
                                    input_prospect = inprospect_green)
  }
  if (sail_version =='4SAIL2'){
    if (is.null(input_prospect)){
      input_prospect <- data.frame('chl' = chl, 'car' = car, 'ant' = ant,
                                   'brown' = brown, 'ewt' = ewt, 'lma' = lma,
                                   'prot' = prot, 'cbc' = cbc,
                                   'n_struct' = n_struct, 'alpha' = alpha)
    }
    # 4SAIL2 requires one of the following combination of input parameters
    # Case #1: optical properties for brown vegetation provided
    if (!is.null(brown_lop)){
      check_brown_lop(brown_lop = brown_lop,
                      lambda = spec_sensor$lambda,
                      input_prospect = input_prospect)
      # Case #2: two sets of input data for prospect
    } else if (is.null(brown_lop)){
      # fraction_brown = 0: green vegetation optics assigned to brown vegetation
      if (fraction_brown==0) {
        brown_lop <- green_lop
      } else {
        if (!dim(input_prospect)[1]==2){
          message('4SAIL2 needs 2 sets of optical props for green & brown veg')
          message('1 set only curently defined. Switch to 4SAIL')
          sail_version <- '4SAIL'
        } else {
          if (utils::packageVersion("prospect")<'2.0.0'){
            inprospect_brown <- prospect::define_Input_PROSPECT(input_prospect[2,])
            brown_lop <- prospect::PROSPECT(SpecPROSPECT = spec_sensor,
                                            Input_PROSPECT = inprospect_brown)
            names(brown_lop) <- c('wvl', 'reflectance', 'transmittance')
          } else {
            inprospect_brown <- prospect::define_input_prospect(input_prospect = input_prospect[2,])
            brown_lop <- prospect::prospect(spec_prospect = spec_sensor,
                                            input_prospect = inprospect_brown)
          }
        }
      }
    }
  }
  return(list('green_lop' = green_lop,
              'brown_lop' = brown_lop))
}


#' @rdname prosail-deprecated
#' @export
adjust_PROSPECT_2_SAIL <- function(SAILversion, Spec_Sensor, Input_PROSPECT,
                                   CHL, CAR, ANT, BROWN, EWT, LMA,
                                   PROT, CBC, N, alpha, fraction_brown,
                                   BrownLOP = NULL){
  .Deprecated("adjust_prospect_to_sail")
  adjust_prospect_to_sail(SAILversion, Spec_Sensor, Input_PROSPECT, CHL, CAR,
                          ANT, BROWN, EWT, LMA, PROT, CBC, N, alpha,
                          fraction_brown, BrownLOP)
}
