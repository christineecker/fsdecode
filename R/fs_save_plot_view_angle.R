
#' fs_save_plot_view_angle.R
#'
#' Saves fsbrain plot with colorbar
#'
#' @param plot fsbrain plot; e.g., fs_plot(fs.overlay)
#' @param write.dir DESCRIPTION.
#' @param filename DESCRIPTION.
#' @param angle.set Image orientation to plot; output of fsbrain::get.view.angle.names()
#' @param colorbar Print colorbar 'yes' or 'no'
#' @param legend.label Colorbar label
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_save_plot_view_angle <- function(plot,
                                    write.dir,
                                    filename,
                                    angle.set = "sd_lateral_lh",
                                    colorbar = "yes",
                                    legend.label = "")
{

  # create output directory --------------------------------------------------

  system(paste0("mkdir ", write.dir))

  # plot image --------------------------------------------------

  brainview_img <- ifelse(colorbar == "yes",
                     paste0(write.dir, "tmp.img.png"),
                     paste0(write.dir, filename))

  fsbrain::vislayout.from.coloredmeshes(
    plot,
    view_angles = angle.set,
    transparency_color = "#FFFFFF",
    rgloptions = list('windowRect' = c(20, 20, 1200, 600)),
    output_img = brainview_img
  )

  ##| print colorbar --------------------------------------------------------
  ##|

  if (colorbar == "yes") {
    fsbrain::coloredmesh.plot.colorbar.separate(
      plot,
      image.plot_extra_options = list(
        'horizontal' = F,
        'legend.lab' = legend.label,
        'legend.width' = 1.1,
        'legend.cex' = 1.8,
        'legend.line' = -3,
        'axis.args' = list(cex.axis = 1.5, family = "Arial")
      ),
      png_options = list(
        'filename' = paste0(write.dir, "tmp.cb.png"),
        width = 430,
        height = 400,
        bg = "#FFFFFF00"
      )
    )

    fsbrain::combine.colorbar.with.brainview.image(
      brainview_img = brainview_img,
      colorbar_img = paste0(write.dir, "tmp.cb.png"),
      output_img = paste0(write.dir, filename),
      horizontal = F,
      offset = "+0+0",
      extend_brainview_img_height_by = 115,
      transparency_color = "#FFFFFF"
    )

    system(paste0("cd ", write.dir, "; rm tmp*"))
  }

  # close images --------------------------------------------------

  while (rgl::rgl.cur() > 0) { rgl::close3d() }

}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- pet.GABAABZ$lh
# hemi <- "lh"
# plot <- fs_plot(fs.overlay, hemi = hemi)
#
# angle.set <- fsbrain::get.view.angle.names(angle_set = "lh")
# angle.set <- "sd_lateral_lh"
#
# fs_save_plot_view_angle(plot,
#                         "output/",
#                         "test.png",
#                         "sd_lateral_lh",
#                         "no",
#                         "test")

