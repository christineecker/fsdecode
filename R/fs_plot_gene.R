#' fs_plot_gene
#'
#' Displays gene expression data for gene
#'
#' @param gene.name DESCRIPTION.
#' @param type Either 'symbol' or 'entrez'
#' @param range description
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' fs_display_gene_expression_map("HTR1A", type = "symbol")
fs_plot_gene <- function(gene.name, type="symbol", range = NULL) {

  if (type == "symbol") {
    gene.index <- which(abagen.genes$symbol == gene.name)
  }

  if (!exists("mRNA")) {
    mRNA <- load_mRNA_expression_maps()
  }

  fs.overlay <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[gene.index,]

  if (is.null(range)) {
    range <- c(min(fs.overlay, na.rm = T) - 0.01,
               max(fs.overlay, na.rm = T) + 0.01)
  }

  cm <- fs.overlay %>%
    fsbrain::vis.data.on.fsaverage(subjects_dir = get_SUBJECTS_DIR(),
                                   vis_subject_id = "fsaverage6",
                                   surface = 'orig',
                                   morph_data_lh = .,
                                   makecmap_options = list('colFn' = matlab::jet.colors,
                                                           'n' = 100,
                                                           'col.na' = "grey",
                                                           'symm'= F,
                                                           'range' = range
                                   ), draw_colorbar = TRUE)

  assign("mRNA", mRNA, envir=.GlobalEnv)
  return(cm)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs_plot_gene("HTR1A", type = "symbol")
