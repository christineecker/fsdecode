#' fs_map2fsavg6
#'
#' Maps freesurfer fsaverage overlay to fsaverage6
#'
#' @param x Vector with FreeSurfer `fsaverage` (fsavg7) overlay data.
#' @param hemi Either "lh" or "rh"
#'
#' @export
#'
#' @return FreeSurfer fsaverage6 (fsavg6) overlay vector.
#' @examples
#' data.fsavg6 <- fs_map2fsavg6(x = data.fsavg7, hemi="lh")
fs_map2fsavg6 <- function(x, hemi="lh") {

  ##| save tmp overlay file --------------------------------------------------------
  ##|

  tmp.filename <- tempfile()
  freesurferformats::write.fs.curv(filepath = tmp.filename, data = fs.overlay)

  ##| specify directories and filenames --------------------------------------------------------
  ##|

  FREESURFER_HOME <-get_FREESURFER_HOME()
  SUBJECTS_DIR <- get_SUBJECTS_DIR()
  trg.filename <- paste0(tmp.filename, ".fsavg6.mgh")

  ##| map to fsavg6 --------------------------------------------------------
  ##|

  script <- paste0(
    "export FREESURFER_HOME=", FREESURFER_HOME, "; cd ", SUBJECTS_DIR, "; export SUBJECTS_DIR=`pwd`; ",
    FREESURFER_HOME, "/bin/mris_apply_reg --src ", tmp.filename,
    " --streg fsaverage/surf/", hemi, ".sphere.reg fsaverage6/surf/", hemi, ".sphere.reg",
    " --trg ", trg.filename)

  system(script)

  ##| read fsavg6 file --------------------------------------------------------
  ##|

  fs.overlay.fsavg6 <- freesurferformats::read.fs.morph(trg.filename)

  ##| remove temp files --------------------------------------------------------
  ##|

  file.remove(tmp.filename, trg.filename)

  return(fs.overlay.fsavg6)
}
