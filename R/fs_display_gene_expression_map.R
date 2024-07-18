#' fs_display_gene_expression_map
#'
#' Displays gene expression data for gene
#'
#' @param gene.name DESCRIPTION.
#' @param type Either 'symbol' or 'entrez'
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' fs_display_gene_expression_map("HTR1A", type = "symbol")
fs_display_gene_expression_map <- function(gene.name, type="symbol") {

  if (type == "symbol") {
    gene.index <- which(abagen.genes$symbol == gene.name)
  }

  if (!exists("mRNA")) {
    mRNA <- load_mRNA_expression_maps()
  }

  cm <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[gene.index,] %>%
    fsbrain::vis.data.on.fsaverage(subjects_dir = get_SUBJECTS_DIR(),
                                   vis_subject_id = "fsaverage6",
                                   surface = 'orig',
                                   morph_data_lh = .,
                                   makecmap_options = list('colFn' = matlab::jet.colors,
                                                           'n' = 100,
                                                           'col.na' = "grey",
                                                           'symm'= F,
                                                           'range' = c(min(., na.rm=T) - 0.01,
                                                                       max(., na.rm=T) + 0.01)
                                   ), draw_colorbar = TRUE)

  assign("mRNA", mRNA, envir=.GlobalEnv)
  return(cm)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs_display_gene_expression_map("HTR1A", type = "symbol")
