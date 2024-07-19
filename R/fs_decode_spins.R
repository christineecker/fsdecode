
#' fs_decode_spins
#'
#' Performs gene expression decoding of FreeSurfer overlay using Alexander-Bloch spin model.
#' Requires that medial wall labels are set to O rather than NA.
#'
#' @param fs.overlay FreeSurfer fsaverage6 overlay for both hemispheres
#' @param data.env Data environment with mRNA data
#' @param n.spins Number of spins
#' @param genes Gene list in symbols
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_decode_spins <- function(fs.overlay,
                            data.env = NULL,
                            n.spins = 1000L,
                            genes = abagen.genes$symbol)
{

  ##| create data environment --------------------------------------------------------
  ##|

  fs.overlay[fsavg6$medial.wall.verts$both] <- 0

  data.env <- NULL
  if(is.null(data.env)){
    data.env <- new.env()
    data.env <- fs_create_decode_spins_data_env("both", genes)
  }

  ##| compute spins --------------------------------------------------------
  ##|

  fs.overlay.spins <- fsnulls::fs_create_spins_bloch(fs.overlay, n.spins)
  fs_plot(fs.overlay.spins[2,])

  ##| compute observed correlations  --------------------------------------------------------
  ##|

  source("R/row_corr_eigen.R")
  writeLines("\nComputing observed correlations between overlay and all genes ...")
  r.obs <- vector()
  r.obs <- as.vector(row_cor_eigen(data.env$mRNA, t(fs.overlay))) #| n.genes

  ##| compute null correlations --------------------------------------------------------
  ##|

  writeLines("\nComputing null correlations between overlay spins and genes ...")
  r.null <- matrix()
  r.null <- row_cor_eigen(data.env$mRNA, fs.overlay.spins) #| n.genes x n.spins

  ##| compute permutation p value --------------------------------------------------------
  ##|

  writeLines("\nComputing permutation p-values ...")
  p.perm <- vector()
  for (i in 1:nrow(data.env$mRNA)) {
    ## compute two-tailed permutation p value
    p.perm[i] = (sum(abs(r.null[i,]) >= abs(r.obs[i])) + 1) / (1 + n.spins)
  }

  ##| compute fdr adjusted p-values --------------------------------------------------------
  ##|

  p.fdr.adj <- p.adjust(p.perm, method = "fdr")

  ##| create data.table with results --------------------------------------------------------
  ##|

  gene.index <- which(abagen.genes$symbol %in% genes)
  genes.decoded <- data.table::data.table("gene.symbol" = genes,
                                          "gene.entrez" = abagen.genes$entrez[gene.index],
                                          "r" = r.obs,
                                          "p.perm" = p.perm,
                                          "p.fdr.adj" = p.fdr.adj)


  return(genes.decoded)
}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- c(pet.GABAABZ$lh, pet.GABAABZ$rh)
# data.env <- NULL
# n.spins <- 100L
# genes <- abagen.genes$symbol[1:1000]
#
# genes.decoded <- fs_decode_spins(fs.overlay,
#                                  data.env = data.env,
#                                  n.spins = n.spins,
#                                  genes = genes)
# str(genes.decoded)


