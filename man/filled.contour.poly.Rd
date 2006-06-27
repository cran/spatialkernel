\name{filled.contour.poly}
\alias{filled.contour.poly}
\title{Level (Contour) Plots with Polygonal Boundary}
\description{
  This is a revised version of the \pkg{base} \R function 
  \code{\link{filled.contour}}. It additionally plots a polygonal
  boundary.
}
\usage{
filled.contour.poly(x = seq(min(poly[,1]), max(poly[,1]), 
  len = nrow(z)),
y = seq(min(poly[,2]), max(poly[,2]), len = ncol(z)),
z, poly,
xlim = range(x, finite=TRUE),
ylim = range(y, finite=TRUE),
zlim = range(z, finite=TRUE),
levels = pretty(zlim, nlevels), nlevels = 10,
color.palette = risk.colors,
col = color.palette(length(levels) - 1),
llevels = levels, labels = NULL, labcex = 0.6,
drawlabel = TRUE, method = "flattest",
vfont = c("sans serif", "plain"),
lcol = par("fg"), lty = par("lty"), lwd = par("lwd"),   
plot.title, plot.axes, key.title, key.axes,
asp = NA, xaxs = "i", yaxs = "i", las = 1, axes = TRUE,
...)
}
\arguments{
  \item{x,y}{locations of grid lines at which the values in \code{z} are
    measured. These must be in ascending order. By default, equally
    spaced values from 0 to 1 are used. If \code{x} is a list, its
    components \code{x$x} and \code{x$y} are used for \code{x} and
    \code{y}, respectively. If the list has component \code{z} this is
    used for \code{z}.}
  \item{z}{a matrix containing the values to be plotted (\code{NA}s are
    allowed). Note that \code{x} can be used instead of \code{z} for
    convenience.}
  \item{xlim}{\code{x} limits for the plot.}
  \item{ylim}{\code{y} limits for the plot.}
  \item{zlim}{\code{z} limits for the plot.}  
  \item{poly}{a matrix containing the \code{x,y}-coordinates of the
    vertices of the polygon boundary.}
  \item{levels}{a set of levels which are used to partition the range of
    \code{z}. Must be strictly increasing (and finite). Areas with
    \code{z} values between consecutive levels are painted with the same
    color.}
  \item{nlevels}{if \code{levels} is not specified, the range of \code{z} is
    divided into approximately this many levels.}
  \item{color.palette}{a color palette function used to assign
    colors in the plot.}
  \item{col}{an explicit set of colors to be used in the plot. This
    argument overrides any palette function specification.}
  \item{llevels}{numeric vector of levels at which to draw contour
    lines, default is the same as \code{levels}.}
  \item{labels}{a vector giving the labels for the contour lines. If
    \code{NULL} then the levels are used as labels.}
  \item{labcex}{\code{cex} for contour labelling.}
  \item{drawlabel}{logical, contour lines are labelled if \code{TRUE}.}
  \item{method}{character string specifying where the labels will be
    located. Possible values are "simple", "edge" and "flattest" (the
    default). See the Details section.}
  \item{vfont}{if a character vector of length 2 is specified, then
    Hershey vector fonts are used for the contour labels. The first
    element of the vector selects a typeface and the second element
    selects a fontindex (see \code{text} for more information).}
  \item{lcol}{color for the lines drawn.}
  \item{lty}{line type for the lines drawn.}
  \item{lwd}{line width for the lines drawn.}
  \item{plot.title}{statement which add title to the main plot.}
  \item{plot.axes}{statement which draws axes on the main plot. This
    overrides the default axes.}
  \item{key.title}{statement which adds title to the plot key.}
  \item{key.axes}{statement which draws axes on the plot key. This
    overrides the default axis.}
  \item{asp}{the \code{y/x} aspect ratio, see
    \code{\link{plot.window}}.}
  \item{xaxs}{the \code{x} axis style. The default is to use internal
    labeling.}
   \item{yaxs}{the \code{y} axis style. The default is to use internal
     labeling.}
   \item{las}{the style of labeling to be used. The default is to use
     horizontal labeling.}
   \item{axes, ...}{additional graphical parameters.}
}
\note{
  By defining \code{z} values as \code{NA} at points outside the
  polygonal boundary, \code{filled.contour.poly} produces a contour plot
  within the polygonal boundary.
}
\seealso{
  \code{\link{filled.contour}}, \code{\link{contour}} and 
  \code{\link{pinpoly}}
}
\keyword{hplot}

