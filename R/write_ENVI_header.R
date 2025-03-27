#' writes ENVI hdr file
#'
#' @param HDR content to be written
#' @param hdr_path Path of the hdr file
#'
#' @return None
#' @importFrom stringr str_count
#' @export

write_ENVI_header <- function(HDR, hdr_path) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && stringr::str_count(x, "\\w+") > 1))
      x <- paste0("{", paste(x, collapse = ","), "}")
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(HDR), h, sep = " = ")), con = hdr_path)
  return(invisible())
}
