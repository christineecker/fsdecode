
#' load_surrogate_pcs.R
#'
#' Reads permuted, smoothed re-sampled gene expression maps into the environment.
#'
#' @param none No arguments.
#'
#' @export
#'
#' @return data environment with list of surrogate maps for N=9 PCs with 1000 surrogates each.
#'
#' @examples
#' data.env <- load_surrogate_pcs()
#'
load_surrogate_pcs <- function() {

  data.env <- new.env()
  fpath <- system.file("extdata", "lh.rh.mRNA.pcs.fwhm5.gaussian.sa1000.rda", package="fsdecode")
  load(fpath, envir = data.env)

  return(data.env)
}

##| data.env <- load_surrogate_pcs()
