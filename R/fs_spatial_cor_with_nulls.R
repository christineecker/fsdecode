

#' fs_spatial_cor_with_nulls
#'
#' Computes correlation between two FreeSurfer fsaverage6 overlays using precomputed null models.
#'
#' @param fs.x FreeSurfer fsaverage6 overlay for both hemispheres
#' @param fs.y FreeSurfer fsaverage6 overlay for both hemispheres
#' @param fs.x.nulls Matrix of null models for fs.x with dimensions `n.nulls` x `n.verts`
#' @param method Correlation method, e.g., "pearson", "spearman", "kendall"
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_spatial_cor_with_nulls <- function(fs.x,
                                      fs.y,
                                      fs.x.nulls,
                                      method = "pearson")
{

  ##| Prepare data --------------------------------------------------------
  ##|

  fs.x[fsavg6$medial.wall.verts$both] <- 0
  fs.y[fsavg6$medial.wall.verts$both] <- 0
  fs.x.nulls[, fsavg6$medial.wall.verts$both] <- 0

  ##| compute observed correlation --------------------------------------------------------
  ##|

  source("R/row_cor_eigen.R")
  writeLines("\nComputing observed correlations between overlays ...")
  r.obs <- cor(fs.x, fs.y, use = "everything", method = method)

  ##| compute null correlations --------------------------------------------------------
  ##|

  writeLines("\nComputing null correlations between x spins and y ...")
  r.null <- as.vector(row_cor_eigen(fs.x.nulls, t(fs.y))) #| n.nulls

  ##| compute permutation p value --------------------------------------------------------
  ##|

  writeLines("\nComputing permutation p-values ...")
  p.perm <- (sum(abs(r.null) >= abs(r.obs)) + 1) / (1 + n.nulls)

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

# filename <- system.file("extdata/lh.rh.5HT1A.knn30000.sa1000.rda", package = "fsnulls")
# load(filename)
# fs.x.nulls <- sa.5HT1A
#
# fs_spatial_cor_with_nulls(fs.x = c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh),
#                           fs.y = c(pet.GABAABZ$lh, pet.GABAABZ$rh),
#                           fs.x.nulls = fs.x.nulls,
#                           method = "pearson")

