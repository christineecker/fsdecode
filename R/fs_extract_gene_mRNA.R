

#' fs_extract_gene_mRNA
#'
#' Extracts mRNA data for a given gene on the fsaverage6 surface.
#'
#' @param gene.id DESCRIPTION.
#' @param type Either `symbol` or `entrez`.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
fs_extract_gene_mRNA <- function(gene.id,
                                 type = "symbol")
{

  ##| Identify gene.id --------------------------------------------------------
  ##|

  if (type == "symbol") {
    gene.index <- which(abagen.genes$symbol == gene.id)
  } else if (type == "entrez") {
    gene.index <- which(abagen.genes$entrez == gene.id)
  } else {
    stop("type must be either 'symbol' or 'entrez'")
  }

  ##| load mRNA data into workspace --------------------------------------------------------
  ##|

  if (!exists("mRNA")) {
    mRNA <- load_mRNA_expression_maps()
  }

  ##| extract mRNA data for gene --------------------------------------------------------
  ##|

  gene.mRNA <- list(
    "lh" = mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[gene.index,],
    "rh" = mRNA$lh.rh.mRNA.fsavg6.fwhm5$rh[gene.index,],
    "both" = c(mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[gene.index,],
               mRNA$lh.rh.mRNA.fsavg6.fwhm5$rh[gene.index,])
  )

  return(gene.mRNA)

}
