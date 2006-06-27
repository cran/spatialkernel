\name{pinpoly}
\alias{pinpoly}
\title{Check if Points are within Polygon}
\description{Check the location of point(s) with respect to a polygon.}
\usage{
  pinpoly(poly, pts)
}
\arguments{
\item{poly}{matrix containing the \code{x,y}-coordinates of the
  vertices of the polygon boundary.}
\item{pts}{matrix of containing the \code{x,y}-coordinates of the
  point locations.}
}
\value{
  An integer vector of indicators for each point in \code{pts},
  \item{-1}{error when number of polygon vertices exceeds 3000;}
  \item{0}{outside the polygon;}
  \item{1}{at the polygon boundary;}
  \item{2}{inside the polygon.}
}
\references{
  This Fortran code comes from Wm Randolph Franklin,
  Electrical, Computer, and Systems Engineering Department,
  Rensselaer Polytechnic Institute, Troy, New York, at website 
  \url{http://www.ecse.rpi.edu/Homepages/wrf}.
}
\note{
  This function is provided here so that
  users do not need to load other packages,
  as it is not available in the \pkg{base} \R packages.
  %%THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE
  
  The return values have been changed from the original ones
  so that the point is inside (including at the boundary) if
  positive.
}
\seealso{\code{\link{phat}} and \code{\link{mcseg.test}}}
\keyword{math}
