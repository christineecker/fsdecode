#' fs_create_decode_spins_data_env.R
#'
#' Creates data enviroment for `fs_decode_spins`.
#'
#' @details Note.
#'
#' @param hemi "lh", "rh", or "both". Default is "both"
#' @param genes Gene list of genes to decode as symbols. Default is all abagen genes.
#'
#' @export
#'
#' @return Data environment with `mRNA`, `permuted.pcs`, and `gene.allocation`.
#'
#' @examples
#' data.env <- fs_create_decode_spins_data_env("lh")
#' data.env <- fs_create_decode_spins_data_env("both", abagen.genes$symbol[1:10])
#'
fs_create_decode_spins_data_env <- function(hemi = "both",
                                            genes = abagen.genes$symbol)
{

  ########################### Load data into workspace ###########################

  writeLines("\nLoading mRNA expression maps ... ")
  mRNA <- load_mRNA_expression_maps()$lh.rh.mRNA.fsavg6.fwhm5

  ############################## Create environment ##############################

  gene.index <- which(rownames(mRNA$lh) %in% genes)

  # select desired genes
  mRNA$lh <- mRNA$lh[gene.index,]
  mRNA$rh <- mRNA$rh[gene.index,]

  data.env <- new.env()
  if (hemi == "lh") {
    writeLines("\nCreating data environment for left hemisphere ... ")
    data.env$mRNA <- mRNA$lh
  } else if  (hemi == "rh") {
    writeLines("\nCreating data environment for right hemisphere ... ")
    data.env$mRNA <- mRNA$rh
  } else if  (hemi == "both") {
    writeLines("\nCreating data environment for both hemispheres ... ")
    data.env$mRNA <- cbind(mRNA$lh, mRNA$rh)
  } else {
    writeLines("\nPlease specify either 'lh', 'rh', or 'both'")
  }

  return(data.env)
}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# data.env <- fs_create_decode_spins_data_env("both")
