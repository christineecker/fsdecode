

#' fs_plot
#'
#' Plots FreeSurfer overlay data on fsaverage surface
#'
#' @param fs.overlay DESCRIPTION.
#' @param fsaverage DESCRIPTION.
#' @param hemi Either "lh", "rh", or "both"
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_plot <- function(fs.overlay,
                    fsaverage = "fsaverage6",
                    hemi = "both")
{

  if (hemi == "both") {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = "fsaverage6",
      surface = 'orig',
      morph_data_both = fs.overlay,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = c(min(fs.overlay, na.rm =
                          T) - 0.01, max(fs.overlay, na.rm = T) + 0.01)
      ),
      draw_colorbar = TRUE
    )
  } else if (hemi == "lh") {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = "fsaverage6",
      surface = 'orig',
      morph_data_lh = fs.overlay,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = c(min(fs.overlay, na.rm =
                          T) - 0.01, max(fs.overlay, na.rm = T) + 0.01)
      ),
      draw_colorbar = TRUE
    )
  } else {
    cm <- fsbrain::vis.data.on.fsaverage(
      subjects_dir = get_SUBJECTS_DIR(),
      vis_subject_id = "fsaverage6",
      surface = 'orig',
      morph_data_rh = fs.overlay,
      makecmap_options = list(
        'colFn' = colorRamps::matlab.like,
        'n' = 100,
        # provide n to color function
        'col.na' = "grey",
        'symm' = F,
        'range' = c(min(fs.overlay, na.rm =
                          T) - 0.01, max(fs.overlay, na.rm = T) + 0.01)
      ),
      draw_colorbar = TRUE
    )
  }

  return(cm)
}

