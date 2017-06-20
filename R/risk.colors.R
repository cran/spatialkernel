#' Color Palette
#' 
#' This color palette is designed to show risk levels with different colours 
#'   from small risk (bright orange) to high risk (dark red).
#' @param n number of colors (>= 1) to be in the palette.
#' @seealso \code{\link{metre}}, \code{\link{colors}}, and \code{\link{palette}}.
#' @examples
#'   ## risk pie with ten levels
#'   pie(rep(1,10), labels = seq(0.1, 1, 0.1), col = risk.colors(10))
#' @keywords color graphs
#' @export
risk.colors <- function(n)
{
    j <- n%/%4
    c(rgb(0.86*(1:j)/j, 0, 0), 
      rainbow(n-2*j, start=1/20,end=1/7),
      hsv(h = 1/6,s = seq(from = 1 - 1/(2 * j), to = 1/(2 * j), length = j), v = 1))[n:1]
}

