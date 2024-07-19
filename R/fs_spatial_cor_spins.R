

#' fs_spatial_cor_spins
#'
#' Computes spatial correlation between two FreeSurfer fsaverage6 overlays with permutation p-value
#' based on Alexander-Bloch spins
#'
#' @param fs.x Freesurfer fsaverage6 overlay
#' @param fs.y Freesurfer fsaverage6 overlay
#' @param n.spins Number of spins
#' @param method Correlation method, e.g., "pearson", "spearman", "kendall"
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_spatial_cor_spins <- function(fs.x,
                                 fs.y,
                                 n.spins = 100L,
                                 method = "pearson")
{

  ##| generate spins --------------------------------------------------------
  ##|

  fs.x[fsavg6$medial.wall.verts$both] <- 0
  fs.y[fsavg6$medial.wall.verts$both] <- 0

  fs.x.spins <- fsnulls::fs_create_spins_bloch(fs.x, n.spins) #| n.spins x n.verts
  fs_plot(fs.x.spins[2,])

  ##| compute observed correlation --------------------------------------------------------
  ##|

  source("R/row_corr_eigen.R")
  writeLines("\nComputing observed correlations between overlays ...")
  r.obs <- cor(fs.x, fs.y, use = "everything", method = method)
  r.obs

  ##| compute null correlations --------------------------------------------------------
  ##|

  source("R/row_corr_eigen.R")
  writeLines("\nComputing null correlations between x spins and y ...")
  r.null <- as.vector(row_cor_eigen(fs.x.spins, t(fs.y))) #| n.spins

  ##| compute permutation p value --------------------------------------------------------
  ##|

  writeLines("\nComputing permutation p-values ...")
  p.perm <- (sum(abs(r.null) >= abs(r.obs)) + 1) / (1 + n.spins)

  ##| create data.table with results --------------------------------------------------------
  ##|

  results <- data.table::data.table(
    "r.obs" = r.obs,
    "p.perm" = p.perm
  )

  return(results)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.x <- c(pet.GABAABZ$lh, pet.GABAABZ$rh)
# fs.y <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
#
# fs_spatial_cor_spins(fs.x,
#                      fs.y,
#                      n.spins = 100L,
#                      method = "pearson")

