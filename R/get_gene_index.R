

#' get_gene_index
#'
#' Returns gene index in `abagen.genes` list based on symbol
#'
#' @param gene Gene symbol
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#'
get_gene_index <- function(gene) {

  index <- which(abagen.genes$symbol == gene)
  print(index)

  return(index)
}
