.spatialkernelOptions <- new.env(FALSE, globalenv())
## 1--gaussian; 2--quadratic(Epanechnikov); 3--quartic; 
#.adaptpara <- list(kernel = 1, PACKAGE="spatialkerenl")
assign(".adaptpara", list(kernel = 1), envir=.spatialkernelOptions)
assign("kernames", c("gaussian", "epanechnikov", "quartic"), envir=.spatialkernelOptions)
assign("ker4names", c(get("kernames", envir=.spatialkernelOptions), "quadratic"), envir=.spatialkernelOptions) ## equal to "ep"

## check .adaptpara in .spatialkernelOptions for existence and validation 
chkernel <- function()
{ 
  adapt <- get(".adaptpara", envir=.spatialkernelOptions)
	if(!is.list(adapt)) {
	  stop(".adaptpara is reserved for spatialkernel internal usage.") 
  }
  adapt
}

#' Select Smoothing Kernel Function
#' 
#' Select a kernel function for kernel regression and kernel smoothing.
#' @param kernel character string giving the smoothing kernel to be used. 
#'     This must be one of \emph{gaussian}, \emph{epanechnikov}, \emph{quadratic},
#'     \emph{quartic}, or \code{NULL}, and may be abbreviated to a unique prefix.
#' @return A character string of the kernel function selected, or the kernel
#'   function currently being used when \code{kernel} is \code{NULL}.
#' @note The default kernel used is \emph{Gaussian}. Unless users want to use a
#'   non-default kernel, there is no need to call \code{setkernel}.
#'   \emph{quadratic} is an alias for \emph{epanechnikov}.
#'   
#'   \code{setkernel} setup kernel function for both kernel regression in the
#'   type-specific probability estimation and the kernel smoothing in the
#'   intensity function estimation.
#' @seealso \code{\link{cvloglk}}, \code{\link{phat}} and \code{\link{lambdahat}}
#' @examples
#' \dontrun{
#'   setkernel("e") ## Select "epanechnikov" kernel
#'   setkernel()    ## show the kernel currrently being used
#' }
#' @keywords distribution smooth
#' @export
setkernel <- function(kernel=NULL)
{
  adapt <- chkernel()
  if(is.null(kernel)) {
    kf <- get("kernames", envir=.spatialkernelOptions)[adapt$kernel]
  } else {
    kernel <- tolower(kernel)
    kernel <- match.arg(kernel, get("ker4names", envir=.spatialkernelOptions))
    adapt$kernel = switch(kernel,
	  gaussian = 1,
	  quadratic = 2,
	  epanechnikov = 2,
	  quartic = 3)
    assign(".adaptpara", adapt, envir=.spatialkernelOptions)
    kf <- kernel
  }
  kf
}
