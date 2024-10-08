
#' fs_decode_with_nulls
#'
#' Performs gene expression decoding of FreeSurfer fsaverage6 overlay based on precomputed null maps
#' Note. requires null matrix with dimensions n.nulls x n.verts
#'
#' @param fs.overlay FreeSurfer fsaverage6 overlay for both hemispheres
#' @param fs.surrogates Matrix of null models with dimensions `n.nulls` x `n.verts`
#' @param data.env Data environment with mRNA data
#' @param genes Gene list in symbols
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_decode_with_nulls <- function(fs.overlay,
                                 fs.nulls = NULL,
                                 data.env = NULL,
                                 genes = abagen.genes$symbol)
{
  ##| specify for testing --------------------------------------------------------
  ##|

  # fs.overlay <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
  # data.env <- NULL
  # genes <- abagen.genes$symbol
  # fs.surrogates <- spins

  ##| create data environment --------------------------------------------------------
  ##|

  if(is.null(data.env)){
    data.env <- new.env()
    data.env <- fs_create_decode_spins_data_env("both", genes) #| $mRNA: n.genes x n.verts
  }

  if(is.null(fs.nulls)){
    writeLines("\nRequires surrogate maps with dimensions n.perm x n.verts ...")
    break
  }

  ##| Set NAs to zero --------------------------------------------------------
  ##|

  fs.overlay[is.na(fs.overlay)] <- 0
  fs.nulls[is.na(fs.nulls)] <- 0
  data.env$mRNA[is.na(data.env$mRNA)] <- 0

  ##| compute observed correlations  --------------------------------------------------------
  ##|

  source("R/row_cor_eigen.R")
  writeLines("\nComputing observed correlations between overlay and all genes ...")
  r.obs <- vector()
  r.obs <- as.vector(row_cor_eigen(data.env$mRNA, t(fs.overlay))) #| 1 x n.genes

  ##| compute null correlations --------------------------------------------------------
  ##|

  writeLines("\nComputing null correlations between overlay spins and genes ...")
  r.null <- matrix()
  r.null <- row_cor_eigen(data.env$mRNA, fs.nulls) #| returns n.genes x n.spins

  ##| compute permutation p value --------------------------------------------------------
  ##|

  n.nulls <- nrow(fs.nulls)
  n.genes <- nrow(data.env$mRNA)

  writeLines("\nComputing permutation p-values ...")
  p.perm <- vector()
  for (i in 1:n.genes) {
    ## compute two-tailed permutation p value
    p.perm[i] = (sum(abs(r.null[i,]) >= abs(r.obs[i])) + 1) / (1 + n.nulls)
  }

  ##| compute maxT corrected p values --------------------------------------------------------
  ##|

  writeLines("\nComputing maxT adjusted p-values ...")
  p.adj <- vector()
  # get maximum absolute correlation across genes within permutations
  maxR <- apply(abs(r.null), 2, function(x) {max(x, na.rm = T)})
  for (i in 1:n.genes) {
    p.adj[i] = (sum(maxR >= abs(r.obs[i])) + 1) / (1 + n.nulls)
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
                                          "p.maxT.adj" = p.adj,
                                          "p.fdr.adj" = p.fdr.adj)

  return(genes.decoded)

}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
#
# nulls <- fsnulls::fs_create_spins_bloch(fs.overlay, 100L)
#
# genes.decoded <- fs_decode_with_nulls (
#   fs.overlay,
#   fs.nulls = nulls,
#   data.env = NULL,
#   genes = abagen.genes$symbol
# )
