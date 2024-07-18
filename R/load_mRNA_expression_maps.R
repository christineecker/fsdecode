
#' load_mRNA_expression_maps.R
#'
#' Loads smoothed vertex-level mRNA expression maps for both hemispheres into data environment.
#'
#' @export
#'
#' @return Data environment with `lh.rh.mRNA.fsavg6.fwhm5` data.
#' @examples
#' mRNA.data <- load_mRNA_expression_maps()
load_mRNA_expression_maps <- function() {

  data.env <- new.env()
  fpath <- system.file("extdata", "lh.rh.mRNA.fsavg6.fwhm5.rda", package="fsdecodev3")
  load(fpath, data.env)

  return(data.env)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# mRNA.data <- load_mRNA_expression_maps()
#
# lh.mRNA <- mRNA.data$lh.rh.mRNA.fsavg6.fwhm5$lh
# rh.mRNA <- mRNA.data$lh.rh.mRNA.fsavg6.fwhm5$rh
#
# gene.index <- 1
# c(lh.data[gene.index,], rh.data[gene.index,]) %>%
# fsbrain::vis.data.on.fsaverage(subjects_dir = get_SUBJECTS_DIR(),
#                               vis_subject_id = "fsaverage6",
#                               surface = 'orig',
#                               morph_data_both = .,
#                               makecmap_options = list('colFn' = colorRamps::matlab.like,
#                                                       'n' = 100, # provide n to color function
#                                                       'col.na' = "grey",
#                                                       'symm'= F,
#                                                       'range' = c(min(., na.rm=T) - 0.01,
#                                                                   max(., na.rm=T) + 0.01)), draw_colorbar = TRUE)
