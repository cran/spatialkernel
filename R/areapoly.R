#' Signed Area of Polygon
#' 
#' Calculate the area of a polygon and its boundary direction. 
#' 
#' @param poly matrix containing the \code{x,y}-coordinates of the
#' vertices of the polygon boundary.
#' @return 
#' \describe{
#'   \item{area}{positive numeric, the area of the polygon.}
#'   \item{sign}{integer, 1 if the polygon orientation is anticlockwise, -1 otherwise.}
#'   \item{poly}{copy of the passed argument \code{poly}.}    
#' }
#' @references Joseph O'Rourke, Computational Geometry in C (2nd Edition),
#'   Cambridge University Press, 2000 edition.
#' @note
#' This function is provided here so that
#' users do not need to load other packages,
#' as it is not available in the base \R packages.
#' @seealso \code{\link{metre}}, \code{\link{risk.colors}}
#' @aliases areapoly
#' @keywords math
#' @export
areapoly <- function(poly) {
## "." in name string may be confused with class method???
    npoly <- nrow(poly)
    asign <- 1
    ans <- .C("area_poly", as.double(poly), as.integer(npoly),
              area=as.double(0), PACKAGE="spatialkernel")$area
    if(ans < 0) {
        asign <- -1
        ans <- -ans
    }
    invisible(list(area=ans, sign=asign, poly=poly))
}
