
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
  fpath <- system.file("extdata/lh.rh.mRNA.fsavg6.fwhm5.rda", package="fsdecode")
  writeLines("\nLoading mRNA maps into data environment ...")
  load(fpath, data.env)

  return(data.env)
}

################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

##| mRNA.data <- load_mRNA_expression_maps()
