#' get_FREESURFER_HOME.R
#'
#' Returns FreeSurfer `FREESURFER_HOME`
#'
#' @param fs.dir Default FreeSurfer installation e.g. `/Applications/freesurfer/`
#'
#' @export
#'
#' @return FreeSurfer `FREESURFER_HOME`
#'
#' @examples
#' get_FREESURFER_HOME(fs.dir = "/Applications/freesurfer/")
#'
get_FREESURFER_HOME <- function(fs.dir = "/Applications/freesurfer/") {

  fs.file.list <- list.files(fs.dir)
  if (any(grepl("ASegStatsLUT.txt", fs.file.list))) {
    fs.version <- ""
  } else {
    fs.version <- list.files(fs.dir)
  }

  FREESURFER_HOME <- paste0(fs.dir, fs.version)

  return(FREESURFER_HOME[1])
}
