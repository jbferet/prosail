#' Deprecated functions in prosail
#'
#' These functions are provided for compatibility with older versions of prosail but are deprecated.
#'
#' @param SAILversion,Spec_Sensor,Input_PROSPECT,lambda See the new functions that replace the old ones
#' @param CHL,CAR,ANT,BROWN,EWT,LMA,PROT,CBC,N,alpha,fraction_brown,BrownLOP See the new functions that replace the old ones
#' @param raster_path,HybridModel,PathOut,SelectedBands,bandname,MaskRaster See the new functions that replace the old ones
#' @param MultiplyingFactor,maxRows,bigRaster,progressBar,filetype See the new functions that replace the old ones
#' @param m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too,t,k,l See the new functions that replace the old ones
#' @param brfMES,brfMOD,xprior,PriorInfoMean,PriorInfoSD,WeightPrior See the new functions that replace the old ones
#' @param InputPROSAIL,SpecPROSPECT,SpecSOIL,SpecATM,BandNames See the new functions that replace the old ones
#' @param HDR,HDRpath,wvl,FWHM,SensorName,LUT,NoiseLevel,NoiseType See the new functions that replace the old ones
#' @param InitialGuess,LowerBound,UpperBound,SpecPROSPECT_Sensor See the new functions that replace the old ones
#' @param SpecATM_Sensor,SpecSOIL_Sensor,TypeLidf,ParmSet,MeritFunction See the new functions that replace the old ones
#' @param TypeDistrib,GaussianDistrib,minval,maxval See the new functions that replace the old ones
#' @param nbSamples,GeomAcq,Codist_LAI See the new functions that replace the old ones
#' @param s2xml,MTD_TL_xml,verbose,atbd,ImPath,InRefl,SRF See the new functions that replace the old ones
#' @param BRF_LUT,AdditiveNoise,MultiplicativeNoise See the new functions that replace the old ones
#' @param InputVar,nbEnsemble,WithReplacement,method,RegressionModels,Refl See the new functions that replace the old ones
#' @param rdot,rsot,tts,skyl,rsdstar,rddstar,PAR_range,abs_dir,abs_hem See the new functions that replace the old ones
#' @param xinit,parms_xinit,Parms2Estimate,Parm2Set,InVar,Parms2Prior See the new functions that replace the old ones
#' @param LIDFa,LIDFb,q,tto,psi,rsoil,diss,Cv,Zeta See the new functions that replace the old ones
#' @param SpectralProps,Path_SensorResponse,SaveSRF See the new functions that replace the old ones
#'
#'
#' @rdname deprecated
#' @name deprecated
#'
NULL

#' @rdname deprecated
#' @export
adjust_PROSPECT_2_SAIL <- function(SAILversion, Spec_Sensor, Input_PROSPECT,
                                   CHL, CAR, ANT, BROWN, EWT, LMA,
                                   PROT, CBC, N, alpha, fraction_brown,
                                   BrownLOP = NULL){
  .Deprecated(old = "adjust_PROSPECT_2_SAIL",
              new = "adjust_prospect_to_sail")
  adjust_prospect_to_sail(SAILversion, Spec_Sensor, Input_PROSPECT, CHL, CAR,
                          ANT, BROWN, EWT, LMA, PROT, CBC, N, alpha,
                          fraction_brown, BrownLOP)
}

#' @rdname deprecated
#' @export
apply_noise_AddMult <- function(BRF_LUT, AdditiveNoise = 0.01,
                                MultiplicativeNoise = 0.02){
  .Deprecated(old = "apply_noise_AddMult",
              new = "apply_noise_addmult")
  apply_noise_addmult(BRF_LUT, AdditiveNoise, MultiplicativeNoise)
}

#' @rdname deprecated
#' @export
Apply_Noise_LUT <- function(LUT, NoiseLevel, NoiseType = 'relative'){
  .Deprecated(old = "Apply_Noise_LUT",
              new = "apply_noise_lut")
  .Deprecated("apply_noise_lut")
  apply_noise_lut(LUT, NoiseLevel, NoiseType)
}

#' @rdname deprecated
#' @export
Apply_prosail_inversion <- function(raster_path, HybridModel, PathOut,
                                    SelectedBands, bandname, MaskRaster = NULL,
                                    MultiplyingFactor = 10000, maxRows = 100,
                                    bigRaster = FALSE, progressBar = TRUE,
                                    filetype = 'GTiff'){
  .Deprecated(old = "Apply_prosail_inversion",
              new = "apply_prosail_inversion")
  options <- list('multiplying_factor' = MultiplyingFactor, 'maxRows' = maxRows,
                  'progressBar' = progressBar, 'filetype' = filetype)

  apply_prosail_inversion(raster_path = raster_path,
                          mask_path = MaskRaster,
                          hybrid_model = HybridModel,
                          output_dir = PathOut,
                          band_names = bandname,
                          selected_bands = SelectedBands,
                          options = options)
}

#' @rdname deprecated
#' @export
applySensorCharacteristics <- function(wvl, InRefl, SRF){
  .Deprecated(old = "applySensorCharacteristics",
              new = "apply_sensor_characteristics")
  apply_sensor_characteristics(wvl, InRefl, SRF)
}

#' @rdname deprecated
#' @export
check_BrownLOP <- function(BrownLOP, lambda, Input_PROSPECT){
  .Deprecated(old = "check_BrownLOP",
              new = "check_brown_lop")
  check_brown_lop(BrownLOP, lambda, Input_PROSPECT)
}

#' @rdname deprecated
#' @export
check_SpectralSampling <- function(SpecPROSPECT, SpecSOIL, SpecATM){
  .Deprecated(old = "check_SpectralSampling",
              new = "check_spectral_sampling")
  check_spectral_sampling(SpecPROSPECT, SpecSOIL, SpecATM)
}

#' @rdname deprecated
#' @export
Compute_albedo  <- function(rsdstar, rddstar, tts, SpecATM_Sensor,
                            PAR_range = c(400, 2400)){
  .Deprecated(old = "Compute_albedo",
              new = "get_albedo")
  get_albedo(rsdstar, rddstar, tts, SpecATM_Sensor, PAR_range)
}

#' @rdname deprecated
#' @export
Compute_BRF <- function(rdot, rsot, tts, SpecATM_Sensor, skyl = NULL){
  .Deprecated(old = "Compute_BRF",
              new = "get_surf_refl")
  get_surf_refl(
    rdot = rdot,
    rsot = rsot,
    tts = tts,
    spec_atm_sensor = SpecATM_Sensor,
    skyl = skyl
  )
}

#' @rdname deprecated
#' @export
Compute_fAPAR  <- function(abs_dir, abs_hem, tts, SpecATM_Sensor,
                           PAR_range = c(400, 700)){
  .Deprecated(old = "Compute_fAPAR",
              new = "get_fapar")
  get_fapar(abs_dir, abs_hem, tts, SpecATM_Sensor, PAR_range)
}

#' @rdname deprecated
#' @export
Compute_SRF <- function(wvl,FWHM, SensorName = 'user_defined'){
  .Deprecated(old = "Compute_SRF",
              new = "get_srf")
  get_srf(wvl = wvl, fwhm = FWHM, sensor_name = SensorName)
}

#' @rdname deprecated
#' @export
ConservativeScattering <- function(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too){
  .Deprecated(old = "ConservativeScattering",
              new = "conservative_scattering")
  conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
}


#' @rdname deprecated
#' @export
CostVal_RMSE_PROSAIL  <- function(brfMES, brfMOD, xprior, PriorInfoMean = NULL,
                                  PriorInfoSD = NULL, WeightPrior = 0.01){
  .Deprecated(old = "CostVal_RMSE_PROSAIL",
              new = "cost_function_rmse_prosail")
  prior_info <- list('mean' = PriorInfoMean,
                     'SD' = PriorInfoSD,
                     'weight_prior' = WeightPrior)
  cost_function_rmse_prosail(brfMES, brfMOD, xprior, prior_info)
}

#' @rdname deprecated
#' @export
Generate_LUT_4SAIL <- function(InputPROSAIL, SpecPROSPECT, SpecSOIL, SpecATM,
                               BandNames = NULL, SAILversion ='4SAIL',
                               BrownLOP = NULL){
  .Deprecated(old = "Generate_LUT_4SAIL",
              new = "generate_lut_4sail")
  generate_lut_4sail(input_prosail = InputPROSAIL,
                     spec_prospect = SpecPROSPECT,
                     spec_soil = SpecSOIL,
                     spec_atm = SpecATM,
                     band_names = BandNames,
                     SAILversion = SAILversion,
                     brown_lop = BrownLOP)
}

#' @rdname deprecated
#' @export
get_HDR_name <- function(ImPath){
  .Deprecated(old = "get_HDR_name",
              new = "get_hdr_name")
  get_hdr_name(ImPath)
}


#' @rdname deprecated
#' @export
get_default_LUT_input <- function(TypeDistrib = NULL,
                                  GaussianDistrib = NULL,
                                  minval = NULL,
                                  maxval = NULL){
  .Deprecated(old = "get_default_LUT_input",
              new = "get_default_lut_input")
  get_default_lut_input(type_distrib = TypeDistrib,
                        gaussian_distrib = GaussianDistrib,
                        minval = minval, maxval = maxval)
}


#' @rdname deprecated
#' @export
get_atbd_LUT_input <- function(nbSamples = 2000, GeomAcq = NULL,
                               Codist_LAI = TRUE){
  .Deprecated(old = "get_atbd_LUT_input",
              new = "get_atbd_lut_input")
  get_atbd_lut_input(nbSamples, GeomAcq, Codist_LAI)
}

#' @rdname deprecated
#' @export
get_InputPROSAIL <- function(atbd = FALSE, GeomAcq = NULL, Codist_LAI = TRUE,
                             minval = NULL, maxval = NULL,
                             TypeDistrib = NULL, GaussianDistrib = NULL,
                             ParmSet = NULL, nbSamples = 2000, verbose = FALSE){
  .Deprecated(old = "get_InputPROSAIL",
              new = "get_input_prosail")
  get_input_prosail(atbd = atbd, geom_acq = GeomAcq,
                    codistribution_lai = Codist_LAI,
                    minval = minval, maxval = maxval,
                    type_distrib = TypeDistrib,
                    gaussian_distrib = GaussianDistrib,
                    parm_set = ParmSet, nb_samples = nbSamples,
                    verbose = verbose)
}

#' @rdname deprecated
#' @export
get_S2geometry <- function(MTD_TL_xml, verbose=FALSE){
  .Deprecated(old = "get_S2geometry",
              new = "get_s2_geometry")
  get_s2_geometry(MTD_TL_xml, verbose)
}

#' @rdname deprecated
#' @export
get_S2geometry_from_SAFE <- function(s2xml){
  .Deprecated(old = "get_S2geometry_from_SAFE",
              new = "get_s2_geometry_from_SAFE")
  get_s2_geometry_from_SAFE(s2xml)
}

#' @rdname deprecated
#' @export
get_S2geometry_from_THEIA <- function(s2xml){
  .Deprecated(old = "get_S2geometry_from_THEIA",
              new = "get_s2_geometry_from_THEIA")
  get_s2_geometry_from_THEIA(s2xml)
}

#' @rdname deprecated
#' @export
GetRadiometry <- function(SensorName = 'user_defined',
                          SpectralProps = NULL,
                          Path_SensorResponse = './',
                          SaveSRF = TRUE){
  .Deprecated(old = "GetRadiometry",
              new = "get_srf_sensor")
  get_srf_sensor(sensor_name = SensorName,
                 wl = SpectralProps$wl,
                 fwhm = SpectralProps$fwhm,
                 srf_path = Path_SensorResponse,
                 save_srf = SaveSRF)
}

#' @rdname deprecated
#' @export
Invert_PROSAIL <- function(brfMES, InitialGuess = NULL, LowerBound, UpperBound,
                           SpecPROSPECT_Sensor, SpecATM_Sensor, SpecSOIL_Sensor,
                           TypeLidf, ParmSet, MeritFunction = "Merit_RMSE_PROSAIL",
                           PriorInfoMean = NULL, PriorInfoSD = NULL,
                           WeightPrior = 0.01) {
  .Deprecated(old = "Invert_PROSAIL",
              new = "invert_prosail")
  prior_info <- list(
    "mean" = PriorInfoMean,
    "SD" = PriorInfoSD,
    "weight_prior" = WeightPrior
  )
  if (is.null(PriorInfoMean) || is.null(PriorInfoSD)){
    prior_info <- NULL
  }

  new_names <- c(
    "CHL" = "chl", "CAR" = "car", "ANT" = "ant", "BROWN" = "brown",
    "EWT" = "ewt",
    "LMA" = "lma", "PROT" = "prot", "CBC" = "cbc", "N" = "n_struct",
    "alpha" = "alpha",
    "LIDFa" = "lidf_a", "LIDFb" = "lidf_b", "lai" = "lai",
    "q" = "hotspot", "tts" = "tts", "tto" = "tto", "psi" = "psi", "psoil" = "psoil"
  )

  names(InitialGuess) <- new_names[names(InitialGuess)]
  names(LowerBound) <- new_names[names(LowerBound)]
  names(UpperBound) <- new_names[names(UpperBound)]
  names(ParmSet) <- new_names[names(ParmSet)]
  if (MeritFunction == "Merit_RMSE_PROSAIL"){
    MeritFunction = "merit_rmse_prosail"
  }
  invert_prosail(
    refl_mes = brfMES,
    initialization = InitialGuess,
    lower_bound = LowerBound,
    upper_bound = UpperBound,
    spec_prospect_sensor = SpecPROSPECT_Sensor,
    spec_atm_sensor = SpecATM_Sensor,
    spec_soil_sensor = SpecSOIL_Sensor,
    type_lidf = TypeLidf,
    parm_set = ParmSet,
    merit_function = MeritFunction,
    prior_info = prior_info
  )
}

#' @rdname deprecated
#' @export
Jfunc1 <- function(k,l,t){
  .Deprecated(old = "Jfunc1",
              new = "jfunc1")
  jfunc1(k,l,t)
}

#' @rdname deprecated
#' @export
Jfunc2 <- function(k,l,t){
  .Deprecated(old = "Jfunc2",
              new = "jfunc2")
  jfunc2(k,l,t)
}

#' @rdname deprecated
#' @export
Jfunc3 <- function(k,l,t){
  .Deprecated(old = "Jfunc3",
              new = "jfunc3")
  jfunc3(k,l,t)
}

#' @rdname deprecated
#' @export
Jfunc4 <- function(m, t){
  .Deprecated(old = "Jfunc4",
              new = "jfunc4")
  jfunc4(m, t)
}

#' @rdname deprecated
#' @export
Merit_RMSE_PROSAIL <- function(xinit, parms_xinit, brfMES, SpecPROSPECT_Sensor,
                               SpecSOIL_Sensor, SpecATM_Sensor, Parms2Estimate,
                               Parm2Set = NULL, ParmSet = NULL, InVar, TypeLidf,
                               PriorInfoMean = NULL, PriorInfoSD = NULL,
                               Parms2Prior = NULL, WeightPrior = 0.01){
  .Deprecated(old = "Merit_RMSE_PROSAIL",
              new = "merit_rmse_prosail")
  prior_info <- list('mean' = PriorInfoMean,
                     'SD' = PriorInfoSD,
                     'weight_prior' = WeightPrior)
  merit_rmse_prosail(xinit, parms_xinit, brfMES, SpecPROSPECT_Sensor,
                     SpecSOIL_Sensor, SpecATM_Sensor, Parms2Estimate,
                     InVar, type_lidf = 2, prior_info = prior_info,
                     parms_to_prior = Parms2Prior)
}

#' @rdname deprecated
#' @export
NonConservativeScattering <- function(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too){
  .Deprecated(old = "NonConservativeScattering",
              new = "non_conservative_scattering")
  non_conservative_scattering(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too)
}

#' @rdname deprecated
#' @export
PrepareSensorSimulation <- function(SpecPROSPECT,SpecSOIL,SpecATM,SRF){
  .Deprecated(old = "PrepareSensorSimulation",
              new = "prepare_sensor_simulation")
  prepare_sensor_simulation(SpecPROSPECT,SpecSOIL,SpecATM,SRF)
}

#' @rdname deprecated
#' @export
PRO4SAIL <- function(Spec_Sensor = NULL, Input_PROSPECT = NULL, N = 1.5,
                     CHL = 40.0, CAR = 8.0, ANT = 0.0, BROWN = 0.0, EWT = 0.01,
                     LMA = NULL, PROT = 0.0, CBC = 0.0, alpha = 40.0,
                     TypeLidf = 2, LIDFa = 60, LIDFb = NULL, lai = 3,
                     q = 0.1, tts = 30, tto = 0, psi = 60, rsoil = NULL,
                     fraction_brown = 0.0, diss = 0.0, Cv = 1, Zeta = 1,
                     SAILversion = '4SAIL', BrownLOP = NULL){
  .Deprecated(old = "PRO4SAIL",
              new = "prosail")
  prosail(spec_sensor = Spec_Sensor, input_prospect = Input_PROSPECT,
          n_struct = N, chl = CHL, car = CAR, ant = ANT, brown = BROWN,
          ewt = EWT, lma = LMA, prot = PROT, cbc = CBC, alpha = alpha,
          type_lidf = TypeLidf, lidf_a = LIDFa, lidf_b = LIDFb, lai = lai,
          hotspot = q, tts = tts, tto = tto, psi = psi, rsoil = rsoil,
          fraction_brown = fraction_brown, diss = diss, cv = Cv, zeta = Zeta,
          SAILversion = SAILversion, brown_lop = BrownLOP)
}

#' @rdname deprecated
#' @export
PROSAIL_Hybrid_Apply <- function(RegressionModels,Refl, progressBar = FALSE){
  .Deprecated(old = "PROSAIL_Hybrid_Apply",
              new = "prosail_hybrid_apply")
  prosail_hybrid_apply(RegressionModels,Refl, progressBar)
}

#' @rdname deprecated
#' @export

PROSAIL_Hybrid_Train <- function(BRF_LUT, InputVar, nbEnsemble = 20,
                                 WithReplacement = FALSE,
                                 method = 'liquidSVM',
                                 verbose = FALSE, progressBar = FALSE){
  .Deprecated(old = "PROSAIL_Hybrid_Train",
              new = "prosail_hybrid_train")
  prosail_hybrid_train(refl_lut = BRF_LUT, input_variables = InputVar,
                       nb_bagg = nbEnsemble, replacement = WithReplacement,
                       method = method, verbose = verbose,
                       progressBar = progressBar)
}

#' @rdname deprecated
#' @export
read_ENVI_header <- function(HDRpath) {
  .Deprecated(old = "read_ENVI_header",
              new = "read_envi_header")
  read_envi_header(hdr_path = HDRpath)
}

#' @rdname deprecated
#' @export
WhichParmPrior <- function(PriorInfoMean, PriorInfoSD) {
  .Deprecated(old = "WhichParmPrior",
              new = "which_parm_prior")
  which_parm_prior(PriorInfoMean, PriorInfoSD)
}

#' @rdname deprecated
#' @export
WhichParameters2Invert <- function(InitialGuess, LowerBound,
                                   UpperBound, ParmSet) {
  .Deprecated(old = "WhichParameters2Invert",
              new = "which_parms_to_invert")
  which_parms_to_invert(InitialGuess, LowerBound, UpperBound, ParmSet)
}

#' @rdname deprecated
#' @export
write_ENVI_header <- function(HDR, HDRpath) {
  .Deprecated(old = "write_ENVI_header",
              new = "write_envi_header")
  write_envi_header(hdr = HDR, hdr_path = HDRpath)
}
