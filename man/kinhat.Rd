\name{kinhat}
\alias{kinhat}
\title{
  Inhomogeneous K-function Estimation
}
\description{
  Estimate the inhomogeneous K function of a non-stationary point pattern.
}
\usage{
  kinhat(pts, lambda, poly, s)
}
\arguments{
  \item{pts}{matrix of the \code{x,y}-coordinates of the point locations.}
  \item{lambda}{intensity function evaluated at the above point locations.}
  \item{poly}{matrix of the \code{x,y}-coordinates of the polygon boundary.}
  \item{s}{vector of distances at which to calculate the K function.}
}
\value{
  A list with components
  \item{k}{values of estimated K at the distances \code{s}.}
  \item{s}{copy of \code{s}.}
}
\details{
  The inhomogeneous K function is a generalization of the usual K
  function defined for a second-order intensity-reweighted stationary
  point process, proposed by Baddeley \emph{et\ al} (2000). 
  
  When the true intensity function is unknown, and is to be estimated
  from the same data as been used to estimate the K function,
  a modified kernel density estimation implemented in \code{\link{lambdahat}}
  with argument \code{gpts=NULL}
  can be used to calculate the estimated intensity at data points.
  See Baddeley \emph{et al} (2000) for details,
  and Diggle, P.J., \emph{et al} (2006) for a cautious note.
}
\note{
  This code is adapted from \pkg{splancs} (Rowlingson and Diggle, 1993)
  fortran code for the estimation of homogeneous K function
  \code{\link[splancs]{khat}}, with edge correction inherited
  for a general polygonal area.
}
\seealso{
  \code{\link[splancs]{khat}}, \code{\link{lambdahat}} 
}
\references{
  \enumerate{
    \item Baddeley, A. J. and Møller, J. and Waagepetersen R.
    (2000) Non and semi-parametric estimation of interaction in 
    inhomogeneous point patterns, \emph{Statistica Neerlandica}, \bold{54}, 
    3, 329--350.
    \item Diggle, P.J., V. G\eqn{\acute{\mathrm{o}}}mez-Rubio,
    P.E. Brown, A.G. Chetwynd and S. Gooding (2006) Second-order
    analysis of inhomogeneous spatial point processes using case-control
    data, \emph{submitted to Biometrics}.
    \item Rowlingson, B. and Diggle, P. (1993) Splancs: spatial point pattern 
    analysis code in S-Plus.  \emph{Computers and Geosciences}, \bold{19}, 
    627--655.
  }
}
\keyword{spatial}
