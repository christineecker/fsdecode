

#' fs_save_plot_both.R
#'
#' Saves left hemisphere fsbrain plot with colorbar
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
fs_save_plot_both <- function(plot, write.dir, filename, legend.label="") {

  # create output directory --------------------------------------------------

  system(paste0("mkdir ", write.dir))

  # plot image --------------------------------------------------

  fsbrain::vislayout.from.coloredmeshes(
    plot,
    view_angles = fsbrain::get.view.angle.names(angle_set = 't4'),
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
      'legend.width' = 1.1, # 1.1
      'legend.cex' = 1.7,
      'legend.line' = -3, #-3
      'legend.mar' = 12, # 12
      'axis.args' = list(cex.axis = 1.5) # 1.5
    ),
    png_options = list(
      'filename' = paste0(write.dir, "tmp.cb.png"),
      width = 400, # 400
      height = 300, # 300
      bg = "#FFFFFF00" # #FFFFFF00
    ),
    trim_png = T
  )

  # combine both --------------------------------------------------

  fsbrain::combine.colorbar.with.brainview.image(
    paste0(write.dir, "tmp.img.png"),
    paste0(write.dir, "tmp.cb.png"),
    paste0(write.dir, filename),
    horizontal = T,
    offset = "+0+310", # 310
    extend_brainview_img_height_by = 5, # 10
    allow_colorbar_shrink = TRUE,
    transparency_color = "#FFFFFF"
  )

  ##| add L and R annotation ===============================================
  ##|
  file <- paste0(write.dir, filename)
  img <- magick::image_read(file)
  img.out <- img %>%
    magick::image_trim(.) %>%
    magick::image_annotate(., "L", size = 30, color = "black", gravity = "West") %>%
    magick::image_annotate(., "R", size = 30, color = "black", gravity = "East")
  magick::image_write(img.out, path = file)

  # close images --------------------------------------------------

  while (rgl::rgl.cur() > 0) { rgl::close3d() }

  ## remove files ------------------------------------------

  system(paste0("cd ", write.dir, "; rm tmp*"))

}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- c(pet.GABAABZ$lh, pet.GABAABZ$rh)
# cm <- fs_plot(fs.overlay)
# fs_save_plot_both(cm, "output/", filename = "test.png")
