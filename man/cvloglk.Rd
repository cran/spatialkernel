\name{cvloglk}
\alias{cvloglk}
\title{Cross-Validated Log-Likelihood Function}
\description{
  Calculate the cross-validated log-likelihood function.
}
\usage{
  cvloglk(pts, marks, t = NULL, h)
}
\details{
  Select a common bandwidth for kernel regression estimation of type-specific
  probabilities of a multivariate Poisson point process with independent 
  component processes of each categorical type by maximizing the cross-validate
  log-likelihood function.
  
  Select a common bandwidth for kernel regression of type-specific 
  probabilities for all time-periods when the argument \code{t} is not 
  \code{NULL}, in which case the data is of a multivariate spatial-temporal 
  point process, with \code{t} the values of associated time-periods.
}
\arguments{
  \item{pts}{matrix containing the \code{x,y}-coordinates of the point 
    locations.}
  \item{marks}{numeric/character vector of the marked labels of the type 
    of each point.}
  \item{t}{numeric vector of the associated time-periods, default
    \code{NULL} for pure spatial data.}
  \item{h}{numeric vector of the kernel smoothing bandwidths at which to
    calculate the cross-validated log-likelihood function.}
}
\value{
  A list with components
  \item{cv}{vector of the values of the cross-validated Log-likelihood
    function.}
  \item{hcv}{numeric value which maximizing the cross-validate log-likelihood
    function}
  \item{...}{copy of the arguments \code{pts, marks, h}.}
}
\references{
  \enumerate{
    \item Diggle, P.J., Zheng, P. and Durr, P. A. (2005)
    Nonparametric estimation of spatial segregation in a multivariate
    point process: bovine tuberculosis in Cornwall, UK. \emph{J. R. 
    Stat. Soc. C}, \bold{54}, 3, 645--658.
  }
}
\seealso{\code{\link{phat}}, \code{\link{mcseg.test}}, and
  \code{\link{mcpat.test}}}
\keyword{smooth}
