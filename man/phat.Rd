\name{phat}
\alias{phat}
\title{Estimate Type-Specific Probabilities}
\description{Estimate the type-specific probabilities for a
  multivariate Poisson point process with independent component
  processes of each type.}
\usage{
phat(gpts, pts, marks, h)
}
\arguments{
  \item{gpts}{matrix containing the \code{x,y}-coordinates of the
    point locations at which type-specific probabilities are estimated.}
  \item{pts}{matrix containing the \code{x,y}-coordinates of the data points.}
  \item{marks}{numeric/character vector of the types of the point in the data.}
  \item{h}{numeric value of the bandwidth used in the kernel regression.}
}

\details{  
  The type-specific probabilities for data \eqn{(x_i, m_i)}, where
  \eqn{x_i} are the spatial point locations and \eqn{m_i} are the 
  categorical mark sequence numbers, \eqn{m_i=1,2,\ldots}, 
  are estimated using the kernel smoothing 
  methodology
  \eqn{\hat p_k(x)=\sum_{i=1}^nw_{ik}(x)I(m_i=k)},
  where \eqn{w_{ik}(x)=w_k(x-x_i)/\sum_{j=1}^n w_k(x-x_j)},
  \eqn{w_k(.)} is the kernel function with bandwidth \eqn{h_k>0},
  \eqn{w_k(x)=w_0(x/h_k)/h_k^2}, and \eqn{w_0(\cdot)} is the
  standardised form of the kernel function. 
  
  The default kernel is the \emph{Gaussian}. Different kernels can be 
  selected by calling \code{\link{setkernel}}.
}
\value{
  A list with components 
  \item{p}{matrix of the type-specific probabilities for all types, with
    the type marks as the matrix row names.}
  \item{...}{copy of the arguments \code{pts, dpts, marks, h}.}
}

\references{
  \enumerate{  
    \item  Diggle, P. J. and Zheng, P. and Durr, P. A. (2005)
    Nonparametric estimation of spatial segregation in a multivariate
    point process: bovine tuberculosis in Cornwall, UK. \emph{J.
    R. Stat. Soc. C}, \bold{54}, 3, 645--658.   
  }
}
\seealso{
  \code{\link{cvloglk}}, \code{\link{mcseg.test}}, and
  \code{\link{setkernel}}.  
}
\keyword{multivariate}
\keyword{spatial}
\keyword{nonparametric}
\keyword{smooth}
\keyword{regression}
