

fs_decode_spins <- function()
{
  # source("R/row_correlation_eigen.R")
  # pet <- c(pet.GABAABZ$lh, pet.GABAABZ$rh)
  # pet[fsavg6$medial.wall.verts$both] <- 0
  #
  # n.spins <- 1000L
  # pet.spins <- fsnulls::fs_create_spins_bloch(pet, n.spins)
  # fs_plot(pet.spins[1,])
  #
  # n.verts <- 81924
  # n.genes <- 1600
  # genes <- matrix(runif(n.genes*n.verts, 0, 1), n.genes, n.verts)


  ##| compute observed correlations  --------------------------------------------------------
  ##|

  writeLines("\nComputing observed correlations between overlay and all genes ...")
  r.obs <- vector()
  r.obs <- as.vector(row_cor_eigen(genes, t(pet))) #| n.genes

  ##| compute null correlations --------------------------------------------------------
  ##|

  source("R/row_correlation_eigen.R")
  writeLines("\nComputing null correlations between overlay spins and genes ...")
  r.null <- matrix()
  r.null <- row_cor_eigen(genes, pet.spins) #| n.genes x n.spins


  ##| compute permutation p value --------------------------------------------------------
  ##|

  writeLines("\nComputing permutation p-values ...")
  p.perm <- vector()
  for (i in 1:n.genes) {
    ## compute two-tailed permutation p value
    p.perm[i] = (sum(abs(r.null[i,]) >= abs(r.obs[i])) + 1) / (1 + n.spins)
  }

  ##| compute fdr adjusted p-values --------------------------------------------------------
  ##|

  p.fdr.adj <- p.adjust(p.perm, method = "fdr")

  ##| create data.table with results --------------------------------------------------------
  ##|



  return(genes.decoded)
}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################
