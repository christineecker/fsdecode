#' get_SUBJECTS_DIR.R
#'
#' Returns FreeSurfer `SUBJECST_DIR`
#'
#' @param fs.dir Default FreeSurfer installation e.g. `/Applications/freesurfer/`
#'
#' @export
#'
#' @return FreeSurfer `SUBJECTS_DIR`
#'
#' @examples
#' get_SUBJECTS_DIR(fs.dir = "/Applications/freesurfer/")
get_SUBJECTS_DIR <- function(fs.dir = "/Applications/freesurfer/") {
  fs.file.list <- list.files(fs.dir)
  if (any(grepl("ASegStatsLUT.txt", fs.file.list))) {
    fs.version <- ""
  } else {
    fs.version <- list.files(fs.dir)
  }
  SUBJECTS_DIR <- paste0(fs.dir, fs.version, "/subjects/")
  return(SUBJECTS_DIR[1])
}
