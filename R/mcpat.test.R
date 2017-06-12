#' Monte Carlo Inference of Temporal Changes in Spatial Segregation
#' An approximate Monte Carlo test of temporal changes in a multivariate
#' spatial-temporal point process.
#' @param pts matrix containing the \code{x,y}-coordinates of the
#'     data point locations.
#' @param marks numeric/character vector of the marked type labels of the
#'     data points.
#' @param t numeric vector of the associated time-periods.
#' @param h numeric vector of the bandwidths at which to calculate the
#'     cross-validated log-likelihood function pooled over times.
#' @param ntest integer with default 100, number of simulations for the Monte 
#'     Carlo test
#' @param proc logical, default \code{TRUE} prints the processing 
#'     messages.
#' @details The spatial-temporal data are denoted as \eqn{(x_i, m_i, t_i)}, 
#'   where \eqn{x_i} are the spatial locations, \eqn{m_i} are the categorical
#'   mark sequence numbers, and \eqn{t_i} are the associated time-periods.
#'   
#'   The null hypothesis is that the type-specific probability surfaces are 
#'   constant over time-periods, \emph{i.e.}, \eqn{p_k(x, t)=p_k(x)}, for any
#'   \eqn{t}, where \eqn{p_k(x, t)} are the type-specific probabilities for
#'   \eqn{k}th category within time-period \eqn{t}.
#'   
#'   Each Monte Carlo simulation is sampled from an approximate \emph{true} 
#'   type-specific probability surfaces --- the estimated one from the data.
#'   Approximately, the simulated data and the original data are samples from
#'   the same probability distribution under the null hypothesis. See Diggle,
#'   P.J. \emph{et al} (2005) for more details.
#' @return A list with components
#' \describe{
#'   \item{pvalue}{\eqn{p}-value of the approximate Monte Carlo test.}
#'   \item{...}{copy of \code{pts, marks, t, h, ntest}.}
#' }
#' @references
#'   \enumerate{
#'     \item Diggle, P. J. and Zheng, P. and Durr, P. A. (2005)
#'     Nonparametric estimation of spatial segregation in a multivariate
#'     point process: bovine tuberculosis in Cornwall, UK. \emph{J.
#'       R. Stat. Soc. C}, \bold{54}, 3, 645--658. 
#'   }
#' @seealso \code{\link{cvloglk}}, \code{\link{phat}}, \code{\link{mcseg.test}}
#' @keywords htest
#' @export
mcpat.test <- function(pts, marks, t, h, ntest=100, proc=TRUE)
{
##Monte Carlo test of change of pattern over time (marks)
  p2k <- function(p, r) { ## find k, sum_{j=1}^{k-1} p_j < r <= sum_{j=1}^k
        k <- 1
        pa <- p[1]
        while(r>pa) {
            k <- k+1
            pa <- pa+p[k]
        }
        k
    }
    
    tpfun <- function(p){
        sum(apply(p, 3, function(x) sum((x-apply(p, c(1,2), mean))^2)))
    }
        
    n <- length(marks)
    mtypes <- sort(unique(marks))
    m <- length(mtypes)
    tt <- sort(unique(t))
    ntt <- length(tt)
    ps <- array(NA, dim=c(n, m, ntt), dimnames = list(NULL, mtypes, NULL)) 
    tp <- NULL
    nh <- length(h)
    if(nh > 1) {
      lcp1 <- cvloglp(pts, marks, t, h)$cv
      hopt <- h[which.max(lcp1)]
    } else hopt <- h
    for(i in 1:ntest){
        if(proc) cat("\rProcessing No.", i, "out of", ntest)
        if(i==1) {
            y1 <- marks
        } else {
            runifn <- runif(n, min=0, max=1)
            for(j in 1:n) {
                y1[j] <- mtypes[p2k(p0[j,], runifn[j])]
            }
        }
        for(j in 1:ntt) {
            ndx <- which(t==tt[j])
            ps[,,j] <- phat(pts, pts[ndx,], y1[ndx], hopt)$p[, mtypes]
        }
        if(i==1) {
            p0 <- apply(ps, 1:2, mean) ## mean of p_j(x, t) over t
        }
        tp <- c(tp, tpfun(ps))
    }
    pv <- (ntest+1-rank(tp)[1])/ntest
    ##pv1 <- (ntest-rank(tp[2:ntest])[1])/(ntest-1)
    invisible(list(pvalue=pv, pts=pts, marks=marks, t=t, h=h, ntest=ntest))
}
