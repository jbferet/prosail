#' make sure of the compatibility of variable names between versions of prosail
#'
#' @param input_prosail list. PROSAIL input variables
#'
#' @return input_prosail numeric. updated input_prosail
#' @export

set_compatibility_parms_to_estimate <- function(parms_to_estimate){
  name_vars <- parms_to_estimate
  old_names <- c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'CBC', 'PROT',
                 'q', 'Cw_rel', 'N', 'Cv', 'Zeta',
                 'LIDFa', 'LIDFb', 'TypeLidf')
  new_names <- c('chl', 'car', 'ant', 'brown', 'ewt', 'lma', 'cbc', 'prot',
                 'hotspot', 'cw_rel', 'n_struct', 'cv', 'zeta',
                 'lidf_a', 'lidf_b', 'type_lidf')
  if (any((old_names %in% parms_to_estimate))){
    update_vars <- which(old_names %in% parms_to_estimate)
    matchval <- match(parms_to_estimate, old_names)
    for (i in seq_len(length(parms_to_estimate)))
      if (!is.na(matchval[i]))
        parms_to_estimate[i] <- new_names[matchval[i]]
  }
  return(parms_to_estimate)
}
