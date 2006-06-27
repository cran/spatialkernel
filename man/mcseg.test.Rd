\name{mcseg.test}
\alias{mcseg.test}
\title{Monte Carlo Test of Spatial Segregation in Multivariate Point Process}
\description{Monte Carlo test of spatial segregation in a multivariate 
  point process by simulating data from random re-labelling of the
  categorical marks.
}
\usage{
  mcseg.test(pts, marks, h, stpts = NULL, ntest = 100, 
    proc = TRUE)
}
\arguments{
  \item{pts}{matrix containing the \code{x,y}-coordinates of the
    data point locations.}
  \item{marks}{numeric/character vector of the marked type labels of the 
    point pattern.}
  \item{h}{numeric vector of the bandwidths at which to calculate
    the cross-validated likelihood function.} 
  \item{stpts}{matrix containing the \code{x,y}-coordinates of the
    locations at which to implement the pointwise segregation test,
    with default \code{NULL} not to do the pointwise segregation test.}
  \item{ntest}{integer with default 100, number of simulations 
    for the Monte Carlo test}
  \item{proc}{logical with default \code{TRUE} to print the processing
    messages.}
}
\details{
  The null hypothesis is that the estimated risk surface is
  spatially constant, \emph{i.e.}, the type-specific probabilities
  are \eqn{p_k(x)=p_k}, for all \eqn{k}, see \code{\link{phat}}. Each
  Monte Carlo simulation is done by relabeling the data categorical
  marks at random
  whilst preserving the observed number of cases of each type. 
 
  The segregation test can also be done pointwise, usually at a fine 
  grid of points, to mark the areas where the estimated type-specific
  probabilities are significantly greater or smaller than the
  spatial average.  
}
\value{
  A list with components
  \item{pvalue}{numeric, \eqn{p}-value of the Monte Carlo test.}
  \item{stpvalue}{matrix, \eqn{p}-values of the test at each point in
    \code{stpts} (if \code{stpts} is not \code{NULL}), with each column corresponds
    to one type}
  \item{...}{copy of the arguments \code{pts, marks, h, stpts, ntest, proc}.}
}
\references{
  \enumerate{
    \item Kelsall, J. E. and Diggle, P. J. (1998) Spatial variation in risk: 
      a nonparametric binary regression approach, \emph{Applied Statistics},
      \bold{47}, 559--573.
    \item Diggle, P. J. and Zheng, P. and Durr, P. A. (2005)
      Nonparametric estimation of spatial segregation in a multivariate
      point process: bovine tuberculosis in Cornwall, UK. \emph{J.
      R. Stat. Soc. C}, \bold{54}, 3, 645--658.   
  }
}
\seealso{\code{\link{cvloglk}} and \code{\link{phat}}}
\keyword{htest}
