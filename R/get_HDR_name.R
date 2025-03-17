#' get hdr name from image file name, assuming it is BIL format
#'
#' @param ImPath path of the image
#'
#' @return corresponding hdr
#' @importFrom tools file_ext file_path_sans_ext
#' @export

get_HDR_name <- function(ImPath) {
  if (tools::file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (tools::file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (tools::file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(tools::file_path_sans_ext(ImPath), ".hdr", sep = "")
  }
  if (!file.exists(ImPathHDR)) {
    message("WARNING : COULD NOT FIND HDR FILE")
    print(ImPathHDR)
    message("Process may stop")
  }
  return(ImPathHDR)
}
