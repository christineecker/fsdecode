

#' close_all_figures
#'
#' Closes all open figures
#'
#' @param  DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
close_all_figures <- function() {
  while (rgl::rgl.cur() > 0) { rgl::rgl.close() }
}
