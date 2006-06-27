\name{setkernel}
\alias{setkernel}
%%\alias{getkernel}
%%\alias{.adaptpara}
\title{Select Smoothing Kernel Function}
\description{
  Select a kernel function for kernel regression and kernel smoothing.
}
\usage{
  setkernel(kernel = NULL)
}
\arguments{
  \item{kernel}{character string giving the smoothing kernel to be used. 
    This must be one of \emph{gaussian}, \emph{epanechnikov}, \emph{quadratic},
    \emph{quartic}, or \code{NULL}, and may be abbreviated to a unique prefix.} 
}
\value{
  A character string of the kernel function selected, or the kernel function
  currently being used when \code{kernel} is \code{NULL}.
}
\note{
  The default kernel used is \emph{Gaussian}.
  Unless users want to use a non-default kernel, there is no need to call 
  \code{setkernel}. \emph{quadratic} is an alias for \emph{epanechnikov}.
  
  \code{setkernel} setup kernel function for both kernel regression in
  the type-specific probability estimation and the kernel smoothing
  in the intensity function estimation.
}
\seealso{\code{\link{cvloglk}}, \code{\link{phat}} and \code{\link{lambdahat}}}
\examples{
\dontrun{
  setkernel("e") ## Select "epanechnikov" kernel
  setkernel()    ## show the kernel currrently being used
}
}
\keyword{distribution}
\keyword{smooth}
