#' get a list of plots defined in a vector file defined as input
#'
#' @param sf_obj list. sf object
#' @param nbdigits numeric. number of digits to define plot names
#'
#' @return list of polygons
#' @importFrom sf st_read
#' @export

get_plot_list_from_sf <- function(sf_obj, nbdigits = 3){
  nbPlots <- length(sf_obj)
  plots <- list()
  for (i in seq_len(nbPlots))
    plots[[i]] <- sf_obj[i]
  names(plots) <- num2char(seq_len(nbPlots), nbdigits = nbdigits)
  return(plots)
}
