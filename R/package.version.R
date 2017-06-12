#' Listing Loaded/Installed Package Versions
#' List version numbers of loaded or installed packages.
#' @details This is a revised version of the \pkg{base} \R function
#'   \code{.package}. It gives both package names and their version numbers.
#' @param all.available logical, if \code{TRUE} return all available packages.
#' @param lib.loc character vector describing the location of \R library trees
#'   to search through, or \code{NULL}. The default value of \code{NULL} 
#'   corresponds to all libraries currently known.
#' @return A list with components
#' \describe{
#'   \item{package}{names of loaded or available packages 
#'     if \code{all.available} is \code{TRUE}.}
#'   \item{version}{associated package versions.}
#' }
#' @seealso \code{\link{.packages}}
#' @keywords utilities
#' @export
package.version <- function(all.available=FALSE, lib.loc=NULL)
{
    ##if(is.null(pkg))
    pkg <- .packages(all.available, lib.loc)
    if(is.null(lib.loc)) libpath <- .Library
    pkgs <- NULL
    vers <- NULL
    for(i in pkg) {
        verline <- NULL
        dnm <- paste(libpath, "/", i, sep="")
        flnm <- paste(dnm, "/DESCRIPTION", sep="")
        if(!file.exists(dnm)) {
            cat("\n", i, "not exists.")
            next
        } else if(!file.exists(flnm)) {
            cat("\nDESCRIPTION for", i, "not found.")
            next
        }
        chlines <- readLines(con=flnm, n=-1, ok=TRUE)
        for(j in chlines) 
            verline <- paste(verline, grep("version:", j, ignore.case=T, value=T))
        version <- gsub("Version: ", "", verline, ignore.case=TRUE)
        version <- gsub(" ", "", version)
        cat("\nVersion of", i, ":", version)
        pkgs <- c(pkgs, i)
        vers <- c(vers, version)
    }
    cat("\n")
    invisible(list(package=pkgs, version=vers))
}
