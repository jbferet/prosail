#' adjust prospect inputs and run prospect in preparation for 4SAIL or 4SAIL2
#'
#' @param SAILversion character. should be '4SAIL' or '4SAIL2'
#' @param Spec_Sensor dataframe. spectral properties required to run PROSPECT
#' @param input_prospect dataframe. includes all prospect input parameters
#' @param CHL numeric. Chlorophyll content (microg.cm-2)
#' @param CAR numeric. Carotenoid content (microg.cm-2)
#' @param ANT numeric. Anthocyain content (microg.cm-2)
#' @param BROWN numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param PROT numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)
#' @param N numeric. Leaf structure parameter
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' @param fraction_brown numeric. fraction of brown vegetation (0-1)
#' @param brown_lop dataframe. brown leaf optical properties, when available
#'
#' (simulation of roughness)
#'
#' @return invisible
#' @export

adjust_PROSPECT_2_SAIL <- function(SAILversion, Spec_Sensor, input_prospect,
                                   CHL, CAR, ANT, BROWN, EWT, LMA, PROT, CBC, N,
                                   alpha, fraction_brown, brown_lop = NULL){

  # for all versions of 4SAIL: get green vegetation
  inprospect_green <- prospect::define_Input_PROSPECT(input_prospect[1,],
                                                      CHL[1], CAR[1], ANT[1],
                                                      BROWN[1], EWT[1], LMA[1],
                                                      PROT[1], CBC[1], N[1],
                                                      alpha[1])
  green_lop <- prospect::PROSPECT(SpecPROSPECT = Spec_Sensor,
                                 Input_PROSPECT = inprospect_green)
  if (SAILversion =='4SAIL2'){
    if (is.null(input_prospect)){
      input_prospect <- data.frame('CHL' = CHL, 'CAR' = CAR, 'ANT' = ANT,
                                   'BROWN' = BROWN, 'EWT' = EWT, 'LMA' = LMA,
                                   'PROT' = PROT, 'CBC' = CBC, 'N' = N,
                                   'alpha' = alpha)
    }
    # 4SAIL2 requires one of the following combination of input parameters
    # Case #1: optical properties for brown vegetation provided
    if (!is.null(brown_lop)){
      check_brown_lop(brown_lop = brown_lop,
                      lambda = Spec_Sensor$lambda,
                      input_prospect = input_prospect)
      # Case #2: two sets of input data for prospect
    } else if (is.null(brown_lop)){
      # fraction_brown = 0: green vegetation optics assigned to brown vegetation
      if (fraction_brown==0) {
        brown_lop <- green_lop
      } else {
        if (!dim(input_prospect)[1]==2){
          message('4SAIL2 needs two sets of optical properties for green and brown vegetation')
          message('Currently one set is defined. will run 4SAIL instead of 4SAIL2')
          SAILversion <- '4SAIL'
        } else {
          inprospect_brown <- prospect::define_Input_PROSPECT(input_prospect[2,])
          brown_lop <- prospect::PROSPECT(SpecPROSPECT = Spec_Sensor,
                                          Input_PROSPECT = inprospect_brown)
        }
      }
    }
  }
  return(list('green_lop' = green_lop,
              'brown_lop' = brown_lop))
}
