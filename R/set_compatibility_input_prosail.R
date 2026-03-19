#' make sure of the compatibility of variable names between versions of prosail
#'
#' @param input_prosail list. PROSAIL input variables
#'
#' @return input_prosail numeric. updated input_prosail
#' @export

set_compatibility_input_prosail <- function(input_prosail){
  name_vars <- names(input_prosail)
  old_names <- c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'CBC', 'PROT',
                 'q', 'Cw_rel', 'N', 'Cv', 'Zeta',
                 'LIDFa', 'LIDFb', 'TypeLidf')
  new_names <- c('chl', 'car', 'ant', 'brown', 'ewt', 'lma', 'cbc', 'prot',
                 'hotspot', 'cw_rel', 'n_struct', 'cv', 'zeta',
                 'lidf_a', 'lidf_b', 'type_lidf')
  if (any((old_names %in% name_vars))){
    update_vars <- which(old_names %in% name_vars)
    for (i in update_vars)
      input_prosail[[new_names[i]]] <- input_prosail[[old_names[i]]]
    message('WARNING: the variable names are deprecated')
    cat(paste0("\033[0;", 31, "m",old_names[update_vars],"\033[0m",", "))
    message('These will be automatically updated to')
    cat(paste0("\033[0;", 32, "m",new_names[update_vars],"\033[0m",", "))
    message('Please update your scripts to avoid errors in next versions')
  }
  return(input_prosail)
}
