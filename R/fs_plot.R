

#' fs_plot
#'
#' Plots FreeSurfer overlay data on fsaverage surface
#'
#' @param fs.overlay DESCRIPTION.
#' @param fsaverage Either `fsaverage6` or `fsaverage`
#' @param hemi Either "lh", "rh", or "both"
#' @param range Either `NULL` or vector with min and max.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_plot <- function(fs.overlay,
                    fsaverage = "fsaverage6",
                    hemi = "both",
                    range = NULL)
{

  if (is.null(range)) {
    range <- c(min(fs.overlay, na.rm = T) - 0.01,
               max(fs.overlay, na.rm = T) + 0.01)
  }

  if (hemi == "both") {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = fsaverage,
      surface = 'orig',
      morph_data_both = fs.overlay,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = range
      ),
      draw_colorbar = TRUE
    )
  } else if (hemi == "lh") {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = fsaverage,
      surface = 'orig',
      morph_data_lh = fs.overlay,
      draw_colorbar = TRUE,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = range)
    )
  } else {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = fsaverage,
      surface = 'orig',
      morph_data_rh = fs.overlay,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = range
      ),
      draw_colorbar = TRUE
    )
  }

  return(cm)
}

