## 1--gaussian; 2--quadratic(Epanechnikov); 3--quartic; 
.adaptpara <- list(kernel = 1, PACKAGE="spatialkerenl")
kernames <- c("gaussian", "epanechnikov", "quartic")
ker4names <- c(kernames, "quadratic") ## equal to "ep"

## check .adaptpara in .GlobalEnv for existence and validation 
chkernel <- function()
{ 
  if(exists(".adaptpara", env=.GlobalEnv)) {
    chk <- FALSE
    adapt <- get(".adaptpara", env=.GlobalEnv)
	if(is.list(adapt)) {
	  if(adapt$PACKAGE != "spatialkerenl") chk = TRUE
	} else chk = TRUE
	if(chk) {
	  stop("\n.adaptpara is reserved for spatialkernel internal usage.\n") 
    }
  } else {
    adapt <- get(".adaptpara", env = getNamespace("spatialkernel"))
  }
  adapt
}
	  
setkernel <- function(kernel=NULL)
{
  adapt <- chkernel()
  if(is.null(kernel)) {
    kf <- kernames[adapt$kernel]
  } else {
    kernel <- tolower(kernel)
    kernel <- match.arg(kernel, ker4names)
    adapt$kernel = switch(kernel,
	  gaussian = 1,
	  quadratic = 2,
	  epanechnikov = 2,
	  quartic = 3)
    assign(".adaptpara", adapt, env=.GlobalEnv)
    kf <- kernel
  }
  kf
}
