#' Check if Points are within Polygon
#' 
#' Check the location of point(s) with respect to a polygon.
#' @param poly matrix containing the \code{x,y}-coordinates of the
#'     vertices of the polygon boundary.
#' @param pts matrix of containing the \code{x,y}-coordinates of the
#'     point locations.
#' @return An integer vector of indicators for each point in \code{pts},
#' \describe{
#'   \item{-1}{error when number of polygon vertices exceeds 3000;}
#'   \item{0}{outside the polygon;}
#'   \item{1}{at the polygon boundary;}
#'   \item{2}{inside the polygon.}
#' }
#' @references This Fortran code comes from Wm Randolph Franklin,
#'   Electrical, Computer, and Systems Engineering Department,
#'   Rensselaer Polytechnic Institute, Troy, New York, at website 
#'   \url{http://www.ecse.rpi.edu/Homepages/wrf}.
#' @note This function is provided here so that users do not need to load other
#'   packages, as it is not available in the \pkg{base} \R packages. THE
#'   VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE
#'   
#'   The return values have been changed from the original ones so that the
#'   point is inside (including at the boundary) if positive.
#' @seealso \code{\link{phat}} and \code{\link{mcseg.test}}
#' @keywords math
#' @export
pinpoly <- function(poly, pts)
{
## -1 outside; 0 on boundary; 1 inside
## number of poly points more than 3000--error
## return values change to -1-error, 0-outside, 1-boundary, 2-inside
  if(nrow(poly)>3000) {
        cat("\nBoundary polygon vertices number exceeds 3000.\n")
        return(-1)
    }
    if(is.matrix(pts)){
        ans<-.Fortran("psnpoly", as.double(pts[,1]), as.double(pts[,2]),
                      as.integer(nrow(pts)), as.double(poly[,1]), as.double(poly[,2]),
		      as.integer(nrow(poly)), 
                      inout=integer(nrow(pts)), PACKAGE="spatialkernel")$inout
    }else{
        ans<-.Fortran("pnpoly", as.double(pts[1]), as.double(pts[2]),
                      as.double(poly[,1]), as.double(poly[,2]), as.integer(nrow(poly)), 
                      inout=as.integer(0), PACKAGE="spatialkernel")$inout
    }
    ans + 1
}
