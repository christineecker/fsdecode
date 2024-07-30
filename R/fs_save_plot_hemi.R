
#' fs_save_plot_hemi.R
#'
#' Saves fsbrain plot with colorbar
#'
#' @param plot fsbrain plot
#' @param write.dir DESCRIPTION.
#' @param filename DESCRIPTION.
#' @param legend.label DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_save_plot_hemi <- function(plot,
                              write.dir,
                              filename,
                              angle.set = "lh",
                              legend.label="")

{

  # create output directory --------------------------------------------------

  system(paste0("mkdir ", write.dir))

  ## adjust offset ------------------------------------------

  #offset <- ifelse (legend.label == '', "+0+30", "+0+70")
  offset <- ifelse (legend.label == '', "+0+20", "+0+20")

  # plot image --------------------------------------------------

  fsbrain::vislayout.from.coloredmeshes(
    plot,
    view_angles = fsbrain::get.view.angle.names(angle_set = angle.set),
    transparency_color = "#FFFFFF",
    rgloptions = list('windowRect' = c(20, 20, 1200, 600)),
    output_img = paste0(write.dir, "tmp.img.png")
  )

  # plot colorbar --------------------------------------------------

  fsbrain::coloredmesh.plot.colorbar.separate(
    plot,
    image.plot_extra_options = list(
      'horizontal' = T,
      'legend.lab' = legend.label,
      'legend.width' = 1.1,
      'legend.cex' = 1.8,
      'legend.line' = -3,
      'axis.args' = list(cex.axis = 1.5,
                         family = "Arial"
                         )
    ),
    png_options = list(
      'filename' = paste0(write.dir, "tmp.cb.png"),
      width = 400,
      height = 300,
      bg = "#FFFFFF00"
    )
  )

  # combine both --------------------------------------------------

  if (angle.set == "lh") {
    height <- 50
  } else {
    height <- 70
  }

  fsbrain::combine.colorbar.with.brainview.image(
    paste0(write.dir, "tmp.img.png"),
    paste0(write.dir, "tmp.cb.png"),
    paste0(write.dir, filename),
    horizontal = T,
    offset = offset,
    extend_brainview_img_height_by = height,
    transparency_color = "#FFFFFF"
  )

  # close images --------------------------------------------------

  while (rgl::rgl.cur() > 0) { rgl::rgl.close() }

  ## remove files ------------------------------------------

  system(paste0("cd ", write.dir, "; rm tmp*"))

}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- pet.GABAABZ$rh
# hemi <- "rh"
# cm <- fs_plot(fs.overlay, hemi = hemi)
# fs_save_plot_hemi(cm, "output/", "test.png", angle.set = hemi, legend.label = "test")
