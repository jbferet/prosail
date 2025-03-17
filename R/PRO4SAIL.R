#' PROSAIL simulation based on a set of combinations of input parameters
#' @param Spec_Sensor list. Includes optical constants required for PROSPECT
#' refractive index, specific absorption coefficients and spectral bands
#' @param Input_PROSPECT  list. PROSPECT input variables
#' @param N numeric. Leaf structure parameter
#' @param CHL numeric. Chlorophyll content (microg.cm-2)
#' @param CAR numeric. Carotenoid content (microg.cm-2)
#' @param ANT numeric. Anthocyain content (microg.cm-2)
#' @param BROWN numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param PROT numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' @param TypeLidf numeric. Type of leaf inclination distribution function
#' @param LIDFa numeric.
#' if TypeLidf ==1, controls the average leaf slope
#' if TypeLidf ==2, controls the average leaf angle
#' @param LIDFb numeric.
#' if TypeLidf ==1, controls the distribution's bimodality
#' if TypeLidf ==2, unused
#' @param lai numeric. Leaf Area Index
#' @param q numeric. Hot Spot parameter
#' @param tts numeric. Sun zeith angle
#' @param tto numeric. Observer zeith angle
#' @param psi numeric. Azimuth Sun / Observer
#' @param rsoil numeric. Soil reflectance
#' @param fraction_brown numeric. Fraction of brown leaf area
#' @param diss numeric. Layer dissociation factor
#' @param Cv numeric. vertical crown cover percentage
#' = % ground area covered with crowns as seen from nadir direction
#' @param Zeta numeric. Tree shape factor
#' = ratio of crown diameter to crown height
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param BrownLOP dataframe. optical properties for brown vegetation
#' - WVL, Reflectance, Transmittance
#' - Set to NULL if use PROSPECT to generate it
#'
#' @return list. rdot,rsot,rddt,rsdt
#' rdot: hemispherical-directional reflectance factor in viewing direction
#' rsot: bi-directional reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rddt: bi-hemispherical reflectance factor
#' @import prospect
#' @export
PRO4SAIL <- function(Spec_Sensor = NULL, Input_PROSPECT = NULL, N = 1.5,
                     CHL = 40.0, CAR = 8.0, ANT = 0.0, BROWN = 0.0, EWT = 0.01,
                     LMA = NULL, PROT = 0.0, CBC = 0.0, alpha = 40.0,
                     TypeLidf = 2, LIDFa = 60, LIDFb = NULL, lai = 3,
                     q = 0.1, tts = 30, tto = 0, psi = 60, rsoil = NULL,
                     fraction_brown = 0.0, diss = 0.0, Cv = 1, Zeta = 1,
                     SAILversion = '4SAIL', BrownLOP = NULL){

  if (is.null(Spec_Sensor)) Spec_Sensor <- prospect::SpecPROSPECT_FullRange
  if (is.null(rsoil)) rsoil <- prosail::SpecSOIL$Dry_Soil
  #	PROSPECT: LEAF OPTICAL PROPERTIES
  LOP <- adjust_PROSPECT_2_SAIL(SAILversion = SAILversion,
                                Spec_Sensor = Spec_Sensor,
                                Input_PROSPECT = Input_PROSPECT,
                                CHL = CHL, CAR = CAR, ANT = ANT, BROWN = BROWN,
                                EWT = EWT, LMA = LMA, PROT = PROT, CBC = CBC,
                                N = N, alpha = alpha, BrownLOP = BrownLOP,
                                fraction_brown = fraction_brown)
  #	SAIL: CANOPY REFLECTANCE
  if (SAILversion == '4SAIL'){
    Ref <- fourSAIL(LeafOptics = LOP$GreenLOP,
                    TypeLidf = TypeLidf, LIDFa = LIDFa, LIDFb = LIDFb,
                    lai = lai, q = q, tts = tts, tto = tto, psi = psi,
                    rsoil = rsoil)
  } else if (SAILversion == '4SAIL2'){
    Ref <- fourSAIL2(leafgreen = LOP$GreenLOP, leafbrown = LOP$BrownLOP,
                     TypeLidf = TypeLidf, LIDFa = LIDFa, LIDFb = LIDFb,
                     lai = lai, q = q, tts = tts, tto = tto, psi = psi,
                     rsoil = rsoil, fraction_brown = fraction_brown,
                     diss = diss, Cv = Cv, Zeta = Zeta)
  }
  return(Ref)
}
