\name{metre}
\alias{metre}
\title{Plot Color Level Metre}
\description{This is a simple function provided here for the convenience
  of users. It adds a key showing how the colors map to the values of 
  an \code{image} plot.}
\usage{
  metre(xl, yb, xr, yt, lab, cols = risk.colors(length(lab) - 1), 
	shift = 0, cex = 1)
}
\arguments{
  \item{xl,yb,xr,yt}{\code{x,y} coordinates of the left-bottom
    and right-top corners of the color level metre bar.}
  \item{lab}{metre level labels in the metre.}
  \item{cols}{associated colours, defaults to use \code{\link{risk.colors}}.}
  \item{shift}{distance to shift the label texts away from the metre bar.}
  \item{cex}{numeric character expansion factor.}
}
%\references{}
\seealso{\code{\link{risk.colors}}}
\keyword{graphs}
\keyword{hplot}
