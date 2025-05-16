#' check if spectral sampling is identical between spec_prospect, spec_soil,
#' spec_atm
#' @param spec_prospect list. Includes optical constants required for PROSPECT
#' @param spec_soil list. Includes either dry soil and wet soil, or a unique
#' soil sample if the psoil parameter is not inverted
#' @param spec_atm list. direct and diffuse radiation for clear conditions
#'
#' @return invisible
#' @export
#'
check_spectral_sampling <- function(spec_prospect, spec_soil, spec_atm){
  l1 <- spec_prospect$lambda
  l2 <- spec_soil$lambda
  l3 <- spec_atm$lambda

  if (!length(l1)==length(l2) | !length(l1)==length(l3) |
      !length(l3)==length(l2)){
    message('matching spectral sampling required: ')
    message('spec_prospect, spec_soil, spec_atm')
    stop()
  } else if (length(unique(l1-l2))>1 | length(unique(l1-l3))>1 |
             length(unique(l3-l2))>1) {
    message('matching spectral sampling required: ')
    message('spec_prospect, spec_soil, spec_atm')
    stop()
  } else if (!unique(l1-l2)==0 | !unique(l1-l3)==0 |
             !unique(l3-l2)==0){
    message('matching spectral sampling required: ')
    message('spec_prospect, spec_soil, spec_atm')
    stop()
  }
  return(invisible())
}
