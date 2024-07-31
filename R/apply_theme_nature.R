
#' apply_theme_nature.R
#'
#' Changes background of ggplot
#'
#' @note To save plot, lead library(ragg), then save using ggsave.
#'
#' @param myggplot DESCRIPTION.
#' @param font DESCRIPTION.
#' @param line.size DESCRIPTION.
#' @param base_size DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
#' ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) %>%
#' theme_fancy(font = "Myriad") +
#' geom_point(size = 1, color = "#2182bc", alpha = 1)
#'
apply_theme_nature <- function(myggplot,
                               base.size = 9,
                               line.size = 0.1,
                               font = "Calibri Light") {

  colfunc <- colorRampPalette(c("white", "grey99", "grey95"))

  g <- grid::rasterGrob(
    colfunc(10),
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    interpolate = TRUE
  )

  pl <- myggplot +
    ggplot2::annotation_custom(
      g,
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
   ggplot2::theme_classic(
      base_family = font,
      base_size = base.size,
      base_line_size = line.size,
      base_rect_size = line.size
    ) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(size = line.size),
                   axis.text = ggplot2::element_text(colour="black"))

  return(pl)
}


#' apply_theme_crisp.R
#'
#' Changes background of ggplot
#'
#' @note To save plot, load library(ragg), then save using ggsave.
#'
#' @param myggplot DESCRIPTION.
#' @param font DESCRIPTION.
#' @param line.size DESCRIPTION.
#' @param base_size DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#'
apply_theme_crisp <- function(myggplot,
                              font = "Arial",
                              line.size = 0.1,
                              base.size = 9) {

  pl <- myggplot +
    ggplot2::theme_linedraw(
      base_family = font,
      base_size = base.size,
      base_line_size = line.size,
      base_rect_size = line.size
    ) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(size = line.size),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

  return(pl)
}

