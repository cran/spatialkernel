#' Plot Color Level Metre
#' 
#' This is a simple function provided here for the convenience of users. It adds
#' a key showing how the colors map to the values of an \code{image} plot.
#' @param xl Left coordinate of the color level metre bar.
#' @param yb Bottom coordinate of the color level metre bar.
#' @param xr Right coordinate of the color level metre bar.
#' @param yt Top coordinate of the color level metre bar.
#' @param lab metre level labels in the metre.
#' @param cols associated colours, defaults to use \code{\link{risk.colors}}.
#' @param shift distance to shift the label texts away from the metre bar.
#' @param cex numeric character expansion factor.
#' @seealso \code{\link{risk.colors}}
#' @keywords graphs hplot
#' @export
metre <- function(xl, yb, xr, yt, lab, cols = risk.colors(length(lab) - 1),
                  shift = 0, cex = 1)
{
##"shift" factor to shift away from the rect legend default=1
    n <- length(lab)-1
    dx <- xr-xl
    dy <- yt-yb
    dxy <- max(dx, dy)/n ##increasing step
    drift <- min(xr-xl, yt-yb)*shift
    if(dx>dy) {
        rect(xl+(1:n-1)*dxy, yb, xl+(1:n)*dxy, yt, col=cols)
        text(xl+(0:n)*dxy, yb-drift, lab, cex=cex, pos=1)
    } else {
        rect(xl, yb+(1:n-1)*dxy, xr, yb+(1:n)*dxy, col=cols)
        text(xr+drift, yb+(0:n)*dxy, lab, cex=cex, pos=4)
    }
	return()
}
