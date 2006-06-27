\name{mcpat.test}
\alias{mcpat.test}
\title{Monte Carlo Inference of Temporal Changes in Spatial Segregation}
\description{
  An approximate Monte Carlo test of temporal changes
  in a multivariate spatial-temporal point process.
}
\usage{
  mcpat.test(pts, marks, t, h, ntest = 100, proc = TRUE)
}
\arguments{
  \item{pts}{matrix containing the \code{x,y}-coordinates of the
    data point locations.}
  \item{marks}{numeric/character vector of the marked type labels of the
    data points.}
  \item{t}{numeric vector of the associated time-periods.}
  \item{h}{numeric vector of the bandwidths at which to calculate the
    cross-validated log-likelihood function pooled over times.}
  \item{ntest}{integer with default 100, number of simulations for the Monte 
    Carlo test}
  \item{proc}{logical, default \code{TRUE} prints the processing 
    messages.}
}
\details{
  The spatial-temporal data are denoted as \eqn{(x_i, m_i, t_i)}, 
  where \eqn{x_i} are the spatial locations, \eqn{m_i} are the 
  categorical mark sequence numbers, and \eqn{t_i} are the associated 
  time-periods.
  
  The null hypothesis is that the type-specific
  probability surfaces are 
  constant over time-periods, \emph{i.e.},
  \eqn{p_k(x, t)=p_k(x)},
  for any \eqn{t}, where \eqn{p_k(x, t)} are the type-specific
  probabilities for \eqn{k}th category within time-period \eqn{t}.

  Each Monte Carlo simulation is sampled from an approximate \emph{true}
  type-specific probability surfaces --- the estimated one from the
  data. Approximately, the simulated data and the original data are
  samples from the same probability distribution under the null
  hypothesis. See Diggle, P.J. \emph{et al} (2005) for more details.
}
\value{
  A list with components
  \item{pvalue}{\eqn{p}-value of the approximate Monte Carlo test.}
  \item{...}{copy of \code{pts, marks, t, h, ntest}.}
}
\references{
  \enumerate{
    \item Diggle, P. J. and Zheng, P. and Durr, P. A. (2005)
    Nonparametric estimation of spatial segregation in a multivariate
    point process: bovine tuberculosis in Cornwall, UK. \emph{J.
      R. Stat. Soc. C}, \bold{54}, 3, 645--658. 
  }
}
\seealso{\code{\link{cvloglk}}, \code{\link{phat}}, \code{\link{mcseg.test}}}
\keyword{htest}
