
#' fs_find_gene_by_pattern
#'
#' FUNCTION_DESCRIPTION
#'
#' @param pattern DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' fs_find_gene_by_pattern(pattern = "HTR")
fs_find_gene_by_pattern <- function(pattern = "HTR") {

  gene.symbols <- abagen.genes$symbol[grepl(pattern, abagen.genes$symbol)]
  print(gene.symbols)

}

# fs_find_gene_by_pattern(pattern = "HTR")
