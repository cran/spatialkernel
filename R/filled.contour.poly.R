#' Level (Contour) Plots with Polygonal Boundary
#' This is a revised version of the \pkg{base} \R function 
#'   \code{\link{filled.contour}}. It additionally plots a polygonal boundary.
#' @param x locations of grid lines at which the values in \code{z} are
#'     measured. These must be in ascending order. By default, equally
#'     spaced values from 0 to 1 are used. If \code{x} is a list, its
#'     components \code{x$x} and \code{x$y} are used for \code{x} and
#'     \code{y}, respectively. If the list has component \code{z} this is
#'     used for \code{z}.
#' @param y See \code{x} arg.
#' @param z a matrix containing the values to be plotted (\code{NA}s are
#'     allowed). Note that \code{x} can be used instead of \code{z} for
#'     convenience.
#' @param xlim \code{x} limits for the plot.
#' @param ylim \code{y} limits for the plot.
#' @param zlim \code{z} limits for the plot.
#' @param poly a matrix containing the \code{x,y}-coordinates of the
#'     vertices of the polygon boundary.
#' @param levels a set of levels which are used to partition the range of
#'     \code{z}. Must be strictly increasing (and finite). Areas with
#'     \code{z} values between consecutive levels are painted with the same
#'     color.
#' @param nlevels if \code{levels} is not specified, the range of \code{z} is
#'     divided into approximately this many levels.
#' @param color.palette a color palette function used to assign
#'     colors in the plot.
#' @param col an explicit set of colors to be used in the plot. This
#'     argument overrides any palette function specification.
#' @param llevels numeric vector of levels at which to draw contour
#'     lines, default is the same as \code{levels}.
#' @param labels a vector giving the labels for the contour lines. If
#'     \code{NULL} then the levels are used as labels.
#' @param labcex \code{cex} for contour labelling.
#' @param drawlabel logical, contour lines are labelled if \code{TRUE}.
#' @param method character string specifying where the labels will be
#'     located. Possible values are "simple", "edge" and "flattest" (the
#'       default). See the Details section.
#' @param vfont if a character vector of length 2 is specified, then
#'     Hershey vector fonts are used for the contour labels. The first
#'     element of the vector selects a typeface and the second element
#'     selects a fontindex (see \code{text} for more information).
#' @param lcol color for the lines drawn.
#' @param lty line type for the lines drawn.
#' @param lwd line width for the lines drawn.
#' @param plot.title statement which add title to the main plot.
#' @param plot.axes statement which draws axes on the main plot. This
#'     overrides the default axes.
#' @param key.title statement which adds title to the plot key.
#' @param key.axes statement which draws axes on the plot key. This
#'     overrides the default axis.
#' @param asp the \code{y/x} aspect ratio, see \code{\link{plot.window}}.
#' @param xaxs the \code{x} axis style. The default is to use internal
#'     labeling.
#' @param yaxs the \code{y} axis style. The default is to use internal
#'     labeling.
#' @param las the style of labeling to be used. The default is to use
#'     horizontal labeling.
#' @param axes Logical. Should axes be drawn? See 
#'   \code{\link[graphics]{plot.default}}.
#' @param ... additional graphical parameters.
#' @note
#'   By defining \code{z} values as \code{NA} at points outside the
#'   polygonal boundary, \code{filled.contour.poly} produces a contour plot
#'   within the polygonal boundary.
#' @seealso \code{\link{filled.contour}}, \code{\link{contour}} and 
#'   \code{\link{pinpoly}}
#' @keywords hplot
#' @name filled.contour.poly 
#' @aliases filled.contour.poly
#' @export
filled.contour.poly <- function (x = seq(min(poly[,1]), max(poly[,1]), len = nrow(z)),
               y = seq(min(poly[,2]), max(poly[,2]), len = ncol(z)), 
               z, poly, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
               zlim = range(z, finite = TRUE), 
               levels = pretty(zlim, nlevels), nlevels = 10,
               color.palette = risk.colors, 
               col = color.palette(length(levels) - 1),
               llevels = levels, labels = NULL, labcex = 0.6,
               drawlabel = TRUE, method = "flattest",
               vfont = c("sans serif", "plain"),
               lcol = par("fg"), lty = par("lty"), lwd = par("lwd"),        
               plot.title, plot.axes, key.title, key.axes, asp = NA, 
               xaxs = "i", yaxs = "i", las = 1, axes = TRUE, ...) 
{
    if(!missing(poly)){
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                } else {
                    z <- x
                    x <- seq(min(poly[,1]), max(poly[,1]), len = nrow(z))
                }
            } else stop("no `z' matrix specified")
        } else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
    } else { ## missing poly
        if (missing(z)) {
            if (!missing(x)&&!missing(y)) {
                poly <- y
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                } else {
                    z <- x
                    x <- seq(min(poly[,1]), max(poly[,1]), len = nrow(z))
                    y <- seq(min(poly[,2]), max(poly[,2]), len = ncol(z))
                }
            } else stop("no `z' and `poly' matrices specified")
        } else if (is.list(x)) {
            poly <- y
            y <- x$y
            x <- x$x
        }          
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing x and y values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
        stop("no proper `z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    #.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
    #                        col = col))
    .filled.contour(as.double(x), as.double(y), z, levels=as.double(levels), 
       col = col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            axis(1)
            axis(2)
        }
    }
    else plot.axes
    contour(x, y, z, nlevels, levels, labels,
            xlim, ylim, zlim,
            labcex, drawlabel, method,
            vfont, axes = FALSE, frame.plot = FALSE,
            lcol,  lty, lwd, add=T)
    polygon(poly)
    box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}
