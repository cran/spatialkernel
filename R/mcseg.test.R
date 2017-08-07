#' Monte Carlo Test of Spatial Segregation in Multivariate Point Process
#' 
#' Monte Carlo test of spatial segregation in a multivariate point process by
#' simulating data from random re-labelling of the categorical marks.
#' @param pts matrix containing the \code{x,y}-coordinates of the data point
#'   locations.
#' @param marks numeric/character vector of the marked type labels of the point
#'   pattern.
#' @param h numeric vector of the bandwidths at which to calculate the
#'   cross-validated likelihood function.
#' @param stpts matrix containing the \code{x,y}-coordinates of the locations at
#'   which to implement the pointwise segregation test, with default \code{NULL}
#'   not to do the pointwise segregation test.
#' @param ntest integer with default 100, number of simulations for the Monte
#'   Carlo test.
#' @param proc logical with default \code{TRUE} to print the processing 
#'   messages.
#' @details The null hypothesis is that the estimated risk surface is spatially
#'   constant, \emph{i.e.}, the type-specific probabilities are
#'   \eqn{p_k(x)=p_k}, for all \eqn{k}, see \code{\link{phat}}. Each Monte Carlo
#'   simulation is done by relabeling the data categorical marks at random 
#'   whilst preserving the observed number of cases of each type.
#'   
#'   The segregation test can also be done pointwise, usually at a fine grid of
#'   points, to mark the areas where the estimated type-specific probabilities
#'   are significantly greater or smaller than the spatial average.
#' @return A list with components
#' \describe{
#'   \item{pvalue}{numeric, \eqn{p}-value of the Monte Carlo test.}
#'   \item{stpvalue}{matrix, \eqn{p}-values of the test at each point in
#'     \code{stpts} (if \code{stpts} is not \code{NULL}), with each column corresponds
#'     to one type}
#'   \item{...}{copy of the arguments \code{pts, marks, h, stpts, ntest, proc}.}
#' }
#' @references
#' \enumerate{ 
#'   \item Kelsall, J. E. and Diggle, P. J. (1998) Spatial variation in risk: a
#'   nonparametric binary regression approach, \emph{Applied Statistics},
#'   \bold{47}, 559--573.
#'   \item Diggle, P. J. and Zheng, P. and Durr, P. A. (2005) Nonparametric
#'   estimation of spatial segregation in a multivariate point process: bovine
#'   tuberculosis in Cornwall, UK. \emph{J. R. Stat. Soc. C}, \bold{54}, 3,
#'   645--658. }
#' @seealso \code{\link{cvloglk}} and \code{\link{phat}}
#' @keywords htest
#' @export
mcseg.test <- function(pts, marks, h, stpts=NULL, ntest=100, proc=TRUE)
{
##spatial variation in the risk surface between different types
##return pointwise tolerance limits (PTLs) for plotting
  adapt <- chkernel()
    ynames <- unique(marks)
    m <- length(ynames)
    ynames0 <- 1:m - 1
    names(ynames0) <- ynames
    y1 <- ynames0[as.character(marks)]
    nh <- length(h)
    stat<-rep(0, ntest)
    n <- length(y1)
    if(!is.null(stpts)) {
        nstpts <- nrow(stpts)
        tct <- matrix(NA, nrow=nstpts*m, ncol=ntest) ## m types????
    }
    c <- matrix(1, nrow=n, ncol=nh)
    alpha <- table(y1)/n
    for(i in 1:ntest){
        if(proc) cat("\rProcessing No.", i, "out of", ntest) 
        if(i==1) y2 <- y1 else y2 <- y1[sample(1:n)]
        if(nh > 1) {
            lcc <- .C("lcn", as.double(pts), as.integer(y2), as.integer(n),
                  as.double(h), as.integer(nh), as.integer(adapt$kernel),
                  as.double(c), lc=double(nh), PACKAGE="spatialkernel")$lc
            ##ophndx <- which(lcc==max(lcc, na.rm=TRUE))
            ##oph <- h[ophndx[1]]
            oph <- h[which.max(lcc)]
        } else oph <- h
        p <- .C("hatpn", as.double(pts), as.integer(n), as.double(pts),
                as.integer(y2), as.integer(n), as.double(oph),
                as.integer(adapt$kernel),
                as.double(c), as.integer(m), p=double(n*m), 
                PACKAGE="spatialkernel")$p
        p <- matrix(p, ncol=m)
        stat[i] <- sum(apply(p, 1, function(x) sum((x-alpha)^2)))
        if(!is.null(stpts)) {
            p <- .C("hatpn", as.double(stpts), as.integer(nstpts), as.double(pts),
                    as.integer(y2), as.integer(n), as.double(oph),
                    as.integer(adapt$kernel),
                    as.double(c), as.integer(m), p=double(nstpts*m), 
                    PACKAGE="spatialkernel")$p
            tct[,i] <- (p-rep(alpha, each=nstpts))
        }
    }
    pv <- (ntest+1-rank(stat)[1])/ntest
    if(!is.null(stpts)) {
        pvtc <- apply(tct, 1, function(x) (ntest+1-rank(x)[1])/ntest)
        pvtc <- matrix(pvtc, ncol=m, dimnames = list(NULL, ynames))
        invisible(list(pvalue=pv, stpvalue=pvtc, pts=pts, marks=marks,
                       h=h, stpts=stpts, ntest=ntest))
    } else {
        invisible(list(pvalue=pv, pts=pts, marks=marks, h=h,
                       ntest=ntest))
    }
}
