#' get_ahba_sample_vertex_neighborhood.R
#'
#' Gets the vertex neighborhood of ahba samples on the fsavg6 or fsavg7 surface.
#'
#' @param fsaverage Either 'fsavg6' or 'fsavg7' (i.e. fsaverage)
#' @param max.distance DESCRIPTION.
#' @param mirror Either "TRUE" or "FALSE". Mirrors left to right and right to left hemisphere samples.
#'
#' @export
#'
#' @return Returns list of samples with corresponding vertex ids within max.distance.
#' @examples
#' sample.vertex.neighborhood <- fsdecode::get_ahba_sample_vertex_neighborhood(fsaverage = "fsavg6", max.distance = 5, mirror = FALSE)
#'
get_ahba_sample_vertex_neighborhood <- function(fsaverage, max.distance = 5, mirror = FALSE) {

  if (fsaverage == "fsavg6") {
    if (isTRUE(mirror)) {
      seed.vertices <- c(abagen.mRNA$fsavg6$vertex_fs_r2l, abagen.mRNA$fsavg6$vertex_fs_l2r)
    } else {
      seed.vertices <- abagen.mRNA$fsavg6$vertex_fs
    }
    surface.lh <- fsavg6$orig$lh
    surface.rh <- fsavg6$orig$rh
  } else if (fsaverage == "fsavg7") {
    if (isTRUE(mirror)) {
      seed.vertices <- c(abagen.mRNA$fsavg6$vertex_fs_r2l, abagen.mRNA$fsavg6$vertex_fs_l2r)
    } else {
      seed.vertices <- abagen.mRNA$fsavg6$vertex_fs
    }
    surface.lh <- fsavg7$orig$lh
    surface.rh <- fsavg7$orig$rh
  } else {
    print("not implemented yet ...")
  }

  n.vertices.hemi <- surface.lh$internal$num_vertices_expected

  writeLines("Computing ahba sample neighborhood ...")
  sample.vertex.neighborhood <- list()
  sample.vertex.neighborhood <- pbmcapply::pbmclapply(seed.vertices, function(i) {
    if (i <= n.vertices.hemi) {
      fsbrain:::geod.vert.neighborhood(surface.lh, i, max_distance = max.distance)$vertices
    } else {
      fsbrain:::geod.vert.neighborhood(surface.rh, i - n.vertices.hemi, max_distance = max.distance)$vertices + n.vertices.hemi
    }
  },
  mc.cores = parallel::detectCores() - 2)

  return(sample.vertex.neighborhood)
}


################################################################################
#                                                                              #
#                                   testing                                    #
#                                                                              #
################################################################################

# sample.vertex.neighborhood <- get_ahba_sample_vertex_neighborhood(fsaverage = "fsavg6", max.distance = 5)
