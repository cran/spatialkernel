#' Estimate Type-Specific Probabilities
#' 
#' Estimate the type-specific probabilities for a multivariate Poisson point
#' process with independent component processes of each type.
#' @param gpts matrix containing the \code{x,y}-coordinates of the
#'     point locations at which type-specific probabilities are estimated.
#' @param pts matrix containing the \code{x,y}-coordinates of the data points.
#' @param marks numeric/character vector of the types of the point in the data.
#' @param h numeric value of the bandwidth used in the kernel regression.
#' @details  The type-specific probabilities for data \eqn{(x_i, m_i)}, where 
#'   \eqn{x_i} are the spatial point locations and \eqn{m_i} are the categorical
#'   mark sequence numbers, \eqn{m_i=1,2,\ldots}, are estimated using the kernel
#'   smoothing methodology \eqn{\hat p_k(x)=\sum_{i=1}^nw_{ik}(x)I(m_i=k)}, 
#'   where \eqn{w_{ik}(x)=w_k(x-x_i)/\sum_{j=1}^n w_k(x-x_j)}, \eqn{w_k(.)} is
#'   the kernel function with bandwidth \eqn{h_k>0}, 
#'   \eqn{w_k(x)=w_0(x/h_k)/h_k^2}, and \eqn{w_0(\cdot)} is the standardised
#'   form of the kernel function.
#'   
#'   The default kernel is the \emph{Gaussian}. Different kernels can be 
#'   selected by calling \code{\link{setkernel}}.
#' @return A list with components 
#' \describe{
#'   \item{p}{matrix of the type-specific probabilities for all types, with
#'     the type marks as the matrix row names.}
#'   \item{...}{copy of the arguments \code{pts, dpts, marks, h}.}
#' }
#' @references
#' \enumerate{ 
#'   \item  Diggle, P. J. and Zheng, P. and Durr, P. A. (2005) 
#'   Nonparametric estimation of spatial segregation in a multivariate point
#'   process: bovine tuberculosis in Cornwall, UK. \emph{J. R. Stat. Soc. C},
#'   \bold{54}, 3, 645--658. 
#' }
#' @seealso
#'   \code{\link{cvloglk}}, \code{\link{mcseg.test}}, and
#'   \code{\link{setkernel}}
#' @keywords multivariate spatial nonparametric smooth regression
#' @export
phat <- function(gpts, pts, marks, h)
{
## {y}={0,1,2,..,m-1} converted inside
## calculate phat at point pts
  adapt <- chkernel()
    ngpts <- length(gpts)/2 ## NO. of points
    n <- length(marks)
    ynames <- names(table(marks))
    m <- length(ynames)
    ynames0 <- 1:m -1 
    names(ynames0) <- ynames
    yy <- ynames0[as.character(marks)]
    c <- rep(1, ngpts)
    ans<-.C("hatpn", as.double(gpts), as.integer(ngpts), as.double(pts),
            as.integer(yy), as.integer(n), as.double(h), 
            as.integer(adapt$kernel),
            as.double(c), as.integer(m), p=double(ngpts*m))$p
    ans <- matrix(ans, ncol=m, dimnames=list(NULL, ynames))
    invisible(list(p=ans, pts=pts, gpts=gpts, marks=marks, h=h))
}
