\name{risk.colors}
\alias{risk.colors}
\title{Color Palette}
\description{
This color palette is designed to show risk levels with different colours 
from small risk (bright orange) to high risk (dark red).
}
\usage{
  risk.colors(n)
}
\arguments{
  \item{n}{number of colors (>= 1) to be in the palette.}
}
\seealso{
  \code{\link{metre}}, \code{\link{colors}}, and \code{\link{palette}}.
}
\examples{
  ## risk pie with ten levels
  pie(rep(1,10), labels = seq(0.1, 1, 0.1), col = risk.colors(10))
}
\keyword{color}
\keyword{graphs}
