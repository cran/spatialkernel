\name{areapoly}
\alias{areapoly}
\title{Signed Area of Polygon}
\description{Calculate the area of a polygon and its 
  boundary direction. 
}
\usage{
areapoly(poly)
}
\arguments{
  \item{poly}{matrix containing the \code{x,y}-coordinates of the
    vertices of the polygon boundary.}
}
\value{
  \item{area}{positive numeric, the area of the polygon.}
  \item{sign}{integer, 1 if the polygon orientation is anticlockwise,
    -1 otherwise.}
  \item{poly}{copy of the passed argument \code{poly}.}    
}
\references{
  Joseph O'Rourke, Computational Geometry in C (2nd Edition),
  Cambridge University Press, 2000 edition.
}
\note{
  This function is provided here so that
  users do not need to load other packages,
  as it is not available in the base \R packages.
}
%%\seealso{\code{\link{metre}}, \code{\link{risk.colors}}}
\keyword{math}
