PROBDIST_CLASS <- "probDist"
#' genProbDist Methods
#' 

#'
#'
#' @export
getSamples <- function(obj, n, ...){
  UseMethod("getSamples", obj)
}

#'
#'
#' @export
fval <- function(obj, x, ...){
  UseMethod("fval", obj)
}

#'
#'
#' @export
Fval <- function(obj, x, ...){
  UseMethod("Fval", obj)
}

#'
#'
#'@export
normalize <- function(obj, ...){
  UseMethod("normalize",obj)
}