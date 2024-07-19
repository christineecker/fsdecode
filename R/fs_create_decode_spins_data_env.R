#' fs_create_decode_spins_data_env.R
#'
#' Creates data enviroment for `fs_decode_spins`.
#'
#' @details Note. vertices in medial wall label are set to zero.
#'
#' @param hemi "lh", "rh", or "both". Default is "both"
#' @param genes Gene list of genes to decode as symbols. Default is all abagen genes.
#'
#' @export
#'
#' @return Data environment with `mRNA` expression data
#'
#' @examples
#' data.env <- fs_create_decode_spins_data_env("lh")
#' data.env <- fs_create_decode_spins_data_env("both", abagen.genes$symbol[1:10])
#'
fs_create_decode_spins_data_env <- function(hemi = "both",
                                            genes = abagen.genes$symbol)
{

  ########################### Load data into workspace ###########################

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
    data.env$mRNA[,fsavg6$medial.wall.verts$lh] <- 0
  } else if  (hemi == "rh") {
    writeLines("\nCreating data environment for right hemisphere ... ")
    data.env$mRNA <- mRNA$rh
    data.env$mRNA[,fsavg6$medial.wall.verts$rh] <- 0
  } else if  (hemi == "both") {
    writeLines("\nConcatenating data across hemispheres ... ")
    data.env$mRNA <- cbind(mRNA$lh, mRNA$rh)
    data.env$mRNA[,fsavg6$medial.wall.verts$both] <- 0
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

# data.env <- fs_create_decode_spins_data_env(hemi = "both",
#                                             genes = abagen.genes$symbol[1:10])
# data.env$mRNA[1:10, 1:10]
