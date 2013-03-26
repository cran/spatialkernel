\name{lambdahat}
\alias{lambdahat}
\title{Kernel Density Estimation of Intensity Function}
\description{Kernel density estimation of the intensity function
  of a two-dimensional point process.}
\usage{
  lambdahat(pts, h, gpts = NULL, poly = NULL, edge = TRUE)
}
\arguments{
  \item{pts}{matrix containing the \code{x,y}-coordinates of the
    data point locations.}
  \item{h}{numeric value of the bandwidth used in the kernel smoothing.}
  \item{gpts}{matrix containing the \code{x,y}-coordinates of point
    locations at which to calculate the intensity function, usually
    a fine grid points within \code{poly}, default \code{NULL} to
    estimate intensity function at data locations.}
  \item{poly}{matrix containing the \code{x,y}-coordinates of the
    vertices of the polygon boundary in an anticlockwise order.}
  \item{edge}{logical, with default \code{TRUE} to do edge-correction.}
}
\details{
  Kernel smoothing methods are widely used to estimate the intensity
  of a spatial point process. One problem which arises
  is the need to handle edge effects. Several methods of 
  edge-correction have been proposed.  
  The adjustment factor proposed in Berman and Diggle (1989)
  is a double integration
  \eqn{\int_AK[(x-x_0)/h]/h^2}{int_AK[(x-x_0)/h]/h^2}, where 
  \eqn{A} is a polygonal area, \eqn{K} is the smoothing kernel 
  and \eqn{h} is the bandwidth used for the smoothing.   
  Zheng, P. \emph{et\ al} (2004) proposed an algorithm 
  for fast calculate of Berman and Diggle's adjustment
  factor. 
  
  When \code{gpts} is \code{NULL}, \code{lambdahat} uses a 
  leave-one-out estimator for the intensity at each of the 
  data points, as been suggested in Baddeley \emph{et al} 
  (2000). This leave-one-out estimate at each of the data points
  then can be used in the inhomogeneous K function estimation
  \code{\link{kinhat}} when the true intensity function is unknown.
  
  The default kernel is the \emph{Gaussian}.
  The kernel function is selected by calling \code{\link{setkernel}}.
}
\value{
  A list with components
  \item{lambda}{numeric vector of the estimated intensity function.}
  \item{...}{copy of the arguments \code{pts, gpts, h, poly, edge}.}
}
\references{
  \enumerate{
    \item M. Berman and P. Diggle (1989) Estimating weighted integrals 
      of the second-order intensity of a spatial point process,
      \emph{J. R. Stat. Soc. B}, \bold{51}, 81--92. 
    \item P. Zheng, P.A. Durr and P.J. Diggle (2004) Edge--correction for Spatial
      Kernel Smoothing --- When Is It Necessary? \emph{Proceedings of the GisVet
      Conference 2004}, University of Guelph, Ontario, Canada, June 2004.
	\item Baddeley, A. J. and Møller, J. and Waagepetersen R.
      (2000) Non and semi-parametric estimation of interaction in 
      inhomogeneous point patterns, \emph{Statistica Neerlandica}, \bold{54}, 
      3, 329--350.
    \item  Laurie, D.P. (1982). Algorithm 584 CUBTRI: Adaptive Cubature
      over a Triangle. \emph{ACM--Trans. Math. Software}, \bold{8}, 210--218.
    \item Jonathan R. Shewchuk, \emph{Triangle, a Two-Dimensional Quality Mesh
      Generator and Delaunay Triangulator} at
      \url{http://www-2.cs.cmu.edu/~quake/triangle.html}.
    \item Alan Murta, \emph{General Polygon Clipper} at
      \url{http://www.cs.man.ac.uk/~toby/alan/software/\#gpc}.
    \item NAG's Numerical Library. \emph{Chapter 11: Quadrature}, 
      NAG's Fortran 90 Library.
      (\url{http://www.nag.co.uk/numeric/fn/manual/html/c11_fn03.html})
  }
}
\note{
  In principle, the \emph{double adaptive} double integration algorithm 
  of Zheng, P. \emph{et\ al} (2004) can be applied to other kernel functions. 
  %Furthermore, the area at the present is enclosed by a 
  %simple polygon which could be generalized into a complex area 
  %with polygonal
  %holes inside. For instance, a large lake lays within the land area of study.

  Other source codes used in the implementation of the double integration 
  algorithm include
  \itemize{
    \item Laurie, D.P. (1982) \emph{adaptive cubature} code in Fortran;
    \item Shewchuk, J.R. \emph{triangulation} code in C;
    \item Alan Murta's \emph{polygon intersection} code in C 
      (\emph{Project: Generic Polygon Clipper}).
  }
}
\seealso{\code{\link{setkernel}}, \code{\link{kinhat}}, \code{\link{density}}}
\keyword{nonparametric}
\keyword{smooth}
\keyword{regression}
\keyword{multivariate}
\keyword{spatial}
