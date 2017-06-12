#' Cross-Validated Log-Likelihood Function
#' Calculate the cross-validated log-likelihood function.
#' 
#' Select a common bandwidth for kernel regression estimation of type-specific
#' probabilities of a multivariate Poisson point process with independent
#' component processes of each categorical type by maximizing the cross-validate
#' log-likelihood function.
#' 
#' Select a common bandwidth for kernel regression of type-specific 
#' probabilities for all time-periods when the argument \code{t} is not 
#' \code{NULL}, in which case the data is of a multivariate spatial-temporal 
#' point process, with \code{t} the values of associated time-periods.
#' @param pts matrix containing the \code{x,y}-coordinates of the point
#'   locations.
#' @param marks numeric/character vector of the marked labels of the type of
#'   each point.
#' @param t numeric vector of the associated time-periods, default \code{NULL}
#'   for pure spatial data.
#' @param h numeric vector of the kernel smoothing bandwidths at which to
#'   calculate the cross-validated log-likelihood function.
#' @return 
#' A list with components
#' \describe{
#'   \item{cv}{vector of the values of the cross-validated Log-likelihood
#'     function.}
#'   \item{hcv}{numeric value which maximizing the cross-validate log-likelihood
#'   function}
#'   \item{...}{copy of the arguments \code{pts, marks, h}.}
#' }
#' @references 
#' \enumerate{
#' \item Diggle, P.J., Zheng, P. and Durr, P. A. (2005)
#' Nonparametric estimation of spatial segregation in a multivariate
#'   point process: bovine tuberculosis in Cornwall, UK. \emph{J. R.
#'   Stat. Soc. C}, \bold{54}, 3, 645--658.
#' }
#' @seealso \code{\link{phat}}, \code{\link{mcseg.test}}, and 
#'   \code{\link{mcpat.test}}
#' @aliases cvloglk
#' @keywords smooth
#' @export
cvloglk <- function(pts, marks, t = NULL, h)
{
  if(is.null(t)) {
    ans <- cvlogl(pts, marks, h)
  }	else ans <- cvloglp(pts, marks, t, h)
  ans$hcv <- h[which.max(ans$cv)]
  ans
}

cvlogl <- function(pts, marks, h)
{
## pts[n*2], y[n].
## h[] can be unequally separated values
    ##if(exists(".adaptpara", env=.GlobalEnv)) {
    ##    .adaptpara <- get(".adaptpara", env=.GlobalEnv) ##, env=.GlobalEnv)
    ##} else .adaptpara <- get(".adaptpara", env = getNamespace("spatialkernel"))
	adapt <- chkernel()
    n <- length(marks)
    nh <- length(h)
    types <- unique(marks)
    mtypes <- 1:length(types) - 1 ## y must from 0 to m-1
    names(mtypes) <- types
    y <- mtypes[marks]
    c <- NULL
    for(i in 1:nh) c <- cbind(c, rep(1, n))
    ans<-.C("lcn", as.double(pts), as.integer(y), as.integer(n), as.double(h), 
        as.integer(nh), as.integer(adapt$kernel), as.double(c),
        lc=double(nh))$lc
    invisible(list(cv=ans, pts=pts, marks=marks, h=h))
}

## pooled cvlogl
cvloglp <- function(pts, marks, t, h)
{
    tt <- sort(unique(t))
    ntt <- length(tt)
    lcp <- rep(0, length(h))
    for(i in 1:ntt) {
        ndx <- which(t==tt[i])
        lcp <- lcp+cvlogl(pts[ndx,], marks[ndx], h)$cv
    }
    invisible(list(cv=lcp, pts=pts, marks=marks, t=t, h=h))
}
