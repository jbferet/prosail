#' PROSAIL simulation based on a set of combinations of input parameters
#' @param spec_sensor list. Includes optical constants required for PROSPECT
#' refractive index, specific absorption coefficients and spectral bands
#' @param input_prospect  list. PROSPECT input variables
#' @param n_struct numeric. Leaf structure parameter
#' @param chl numeric. chlorophyll content (microg.cm-2)
#' @param car numeric. carotenoid content (microg.cm-2)
#' @param ant numeric. anthocyain content (microg.cm-2)
#' @param brown numeric. brown pigment content (Arbitrary units)
#' @param ewt numeric. Equivalent Water Thickness (g.cm-2)
#' @param lma numeric. Leaf Mass per Area (g.cm-2)
#' @param prot numeric. protein content  (g.cm-2)
#' @param cbc numeric. NonProtCarbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' @param type_lidf numeric. Type of leaf inclination distribution function
#' @param lidf_a numeric.
#' if type_lidf ==1, controls the average leaf slope
#' if type_lidf ==2, controls the average leaf angle
#' @param lidf_b numeric.
#' if type_lidf ==1, controls the distribution's bimodality
#' if type_lidf ==2, unused
#' @param lai numeric. Leaf Area Index
#' @param hotspot numeric. Hot Spot parameter
#' @param tts numeric. Sun zeith angle
#' @param tto numeric. Observer zeith angle
#' @param psi numeric. Azimuth Sun / Observer
#' @param rsoil numeric. Soil reflectance
#' @param fraction_brown numeric. Fraction of brown leaf area
#' @param diss numeric. Layer dissociation factor
#' @param cv numeric. vertical crown cover percentage
#' = % ground area covered with crowns as seen from nadir direction
#' @param zeta numeric. Tree shape factor
#' = ratio of crown diameter to crown height
#' @param SAILversion character. choose between 4SAIL and 4SAIL2
#' @param brown_lop dataframe. optical properties for brown vegetation
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
prosail <- function(spec_sensor = NULL, input_prospect = NULL, n_struct = 1.5,
                    chl = 40.0, car = 8.0, ant = 0.0, brown = 0.0, ewt = 0.01,
                    lma = NULL, prot = 0.0, cbc = 0.0, alpha = 40.0,
                    type_lidf = 2, lidf_a = 60, lidf_b = NULL, lai = 3,
                    hotspot = 0.1, tts = 30, tto = 0, psi = 60, rsoil = NULL,
                    fraction_brown = 0.0, diss = 0.0, cv = 1, zeta = 1,
                    SAILversion = '4SAIL', brown_lop = NULL){

  if (is.null(spec_sensor))
    spec_sensor <- prosail::spec_prospect_full_range
  if (is.null(rsoil))
    rsoil <- prosail::spec_soil_atbd_v2$soil_01
  #	PROSPECT: LEAF OPTICAL PROPERTIES
  lop <- adjust_prospect_to_sail(sail_version = SAILversion,
                                 spec_sensor = spec_sensor,
                                 input_prospect = input_prospect,
                                 chl = chl, car = car, ant = ant, brown = brown,
                                 ewt = ewt, lma = lma, prot = prot, cbc = cbc,
                                 n_struct = n_struct, alpha = alpha,
                                 brown_lop = brown_lop,
                                 fraction_brown = fraction_brown)
  #	SAIL: CANOPY REFLECTANCE
  if (SAILversion == '4SAIL'){
    refl <- fourSAIL(lop = lop$green_lop,
                     type_lidf = type_lidf, lidf_a = lidf_a, lidf_b = lidf_b,
                     lai = lai, hotspot = hotspot, tts = tts, tto = tto, psi = psi,
                     rsoil = rsoil)
  } else if (SAILversion == '4SAIL2'){
    refl <- fourSAIL2(leaf_green = lop$green_lop, leaf_brown = lop$brown_lop,
                      type_lidf = type_lidf, lidf_a = lidf_a, lidf_b = lidf_b,
                      lai = lai, hotspot = hotspot, tts = tts, tto = tto, psi = psi,
                      rsoil = rsoil, fraction_brown = fraction_brown,
                      diss = diss, cv = cv, zeta = zeta)
  }
  return(refl)
}


#' @rdname prosail-deprecated
#' @export
PRO4SAIL <- function(Spec_Sensor = NULL, Input_PROSPECT = NULL, N = 1.5,
                     CHL = 40.0, CAR = 8.0, ANT = 0.0, BROWN = 0.0, EWT = 0.01,
                     LMA = NULL, PROT = 0.0, CBC = 0.0, alpha = 40.0,
                     TypeLidf = 2, LIDFa = 60, LIDFb = NULL, lai = 3,
                     q = 0.1, tts = 30, tto = 0, psi = 60, rsoil = NULL,
                     fraction_brown = 0.0, diss = 0.0, Cv = 1, Zeta = 1,
                     SAILversion = '4SAIL', BrownLOP = NULL){
  .Deprecated("prosail")
  prosail(Spec_Sensor, Input_PROSPECT, N, CHL, CAR, ANT, BROWN, EWT, LMA, PROT,
          CBC, alpha, TypeLidf, LIDFa, LIDFb, lai, q, tts, tto, psi, rsoil,
          fraction_brown, diss, Cv, Zeta, SAILversion, BrownLOP)
}

