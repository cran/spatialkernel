\name{spatialkernel-package}
\alias{spatialkernel-package}
%\alias{spatialkernel}
\docType{package}
\title{
  The Spatialkernel Package
}
\description{
  An \R package for spatial point process analysis.
}
\details{
This package contains functions for spatial point process
analysis using kernel smoothing methods. This package has
been written to be compatible with the \pkg{splancs} package
which is available on \cite{CRAN} (The Comprehensive R Archive Network).

For a complete list of functions with individual help pages,
use \code{library(help = \ "spatialkernel")}.
}
\author{
  Pingping Zheng and Peter Diggle
}
\section{Maintainer}{
  Pingping Zheng \email{pingping.zheng@lancaster.ac.uk}
}
\references{
  \enumerate{
    \item P. Zheng, P.A. Durr and P.J. Diggle (2004) Edge-correction for Spatial
    Kernel Smoothing --- When Is It Necessary? \emph{Proceedings of the GisVet
      Conference 2004}, University of Guelph, Ontario, Canada, June 2004.
    \item Diggle, P.J., Zheng, P. and Durr, P. A. (2005)
    Nonparametric estimation of spatial segregation in a multivariate
    point process: bovine tuberculosis in Cornwall, UK. \emph{J. R. 
      Stat. Soc. C}, \bold{54}, 3, 645--658.
  }
}
\examples{
## An example of spatial segregation analysis
\dontrun{
## source in Lansing Woods tree data within a polygon boundary
data(lansing)
data(polyb)
## select data points within polygon
ndx <- which(pinpoly(polyb, as.matrix(lansing[c("x", "y")])) > 0)
pts <- as.matrix(lansing[c("x", "y")])[ndx,]
marks <- lansing[["marks"]][ndx]
## select bandwidth
h <- seq(0.02, 0.1, length=101)
cv <- cvloglk(pts, marks, h=h)$cv
hcv <- h[which.max(cv)]
plot(h, cv, type="l")
## estimate type-specific probabilities and do segregation tests
## by one integrated function
sp <- spseg(pts, marks, hcv, opt=3, ntest=1000, poly=polyb)
## plot estimated type-specific probability surfaces
plotphat(sp)
## additional with pointwise significance contour lines
plotmc(sp, quan=c(0.025, 0.975))
## p-value of the Monte Carlo segregation test
cat("\np-value of the Monte Carlo segregation test", sp$pvalue)

##estimate intensity function at grid point for presentation
##with bandwidth hcv
gridxy <- as.matrix(expand.grid(x=seq(0, 1, length=101), y=seq(0, 1, length=101)))
ndx <- which(pinpoly(polyb, gridxy) > 0) ##inside point index
lam <- matrix(NA, ncol=101, nrow=101)
lam[ndx] <- lambdahat(pts, hcv, gpts = gridxy[ndx,], poly =
    polyb)$lambda
brks <- pretty(range(lam, na.rm=TRUE), n=12)
plot(0, 0, xlim=0:1, ylim=0:1, xlab="x", ylab="y", type="n")
image(x=seq(0, 1, length=101), y=seq(0, 1, length=101),
    z=lam, add=TRUE, breaks=brks, col=risk.colors(length(brks)-1))
polygon(polyb)
metre(0, 0.01, 0.05, 0.51, lab=brks, col=risk.colors(length(brks)-1), cex=1)

## An example of inhomogeneous intensity function and K function
## estimated with the same data
s <- seq(0, 0.06, length=101)
lam <- lambdahat(pts, hcv, poly=polyb)$lambda
kin <- kinhat(pts, lam, polyb, s)
plot(kin$s, kin$k-pi*(kin$s)^2, xlab="s", ylab="k-pi*s^2", type="l")
}
}
\note{  
  For the convience of the user, we present here examples which show
  how to use some of the functions in the package.
}
\seealso{\code{\link{cvloglk}}, \code{\link{phat}},
  \code{\link{mcseg.test}}, \code{\link{plotphat}},
  \code{\link{plotmc}}, \code{\link{pinpoly}},
  \code{\link{risk.colors}}, \code{\link{metre}}}
\keyword{ package }
