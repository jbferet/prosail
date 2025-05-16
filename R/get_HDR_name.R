#' get hdr name from image file name, assuming it is BIL format
#'
#' @param image_path path of the image
#'
#' @return corresponding hdr
#' @importFrom tools file_ext file_path_sans_ext
#' @export

get_hdr_name <- function(image_path) {
  if (tools::file_ext(image_path) == "") {
    hdr_path <- paste(image_path, ".hdr", sep = "")
  } else if (tools::file_ext(image_path) == "bil") {
    hdr_path <- gsub(".bil", ".hdr", image_path)
  } else if (tools::file_ext(image_path) == "zip") {
    hdr_path <- gsub(".zip", ".hdr", image_path)
  } else {
    hdr_path <- paste(tools::file_path_sans_ext(image_path), ".hdr", sep = "")
  }
  if (!file.exists(hdr_path)) {
    message("WARNING : COULD NOT FIND HDR FILE")
    print(hdr_path)
    message("Process may stop")
  }
  return(hdr_path)
}
