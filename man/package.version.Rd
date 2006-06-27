\name{package.version}
\alias{package.version}
\title{Listing Loaded/Installed Package Versions}
\description{List version numbers of loaded or installed packages.}
\usage{
  package.version(all.available = FALSE, lib.loc = NULL)
}
\details{
  This is a revised version of the \pkg{base} \R function \code{.package}.
  It gives both package names and their version numbers.
}
\arguments{
  \item{all.available}{logical, if \code{TRUE} return all available packages.}
  \item{lib.loc}{character vector describing the location of \R library 
      trees to search through, or \code{NULL}. The default value of \code{NULL}
      corresponds to all libraries currently known.}
}
\value{
  A list with components
  \item{package}{names of loaded or available packages 
    if \code{all.available} is \code{TRUE}.}
  \item{version}{associated package versions.}
}
\seealso{\code{\link{.packages}}}
\keyword{utilities}
