
#' Cicular geometric median
#'
#' Compute geometric circular median based on the simplified circular distance
#' 
#' @param data Circular dataframe or matrix
#' @param tol geometric median tolerance
#' @param maxiter maximum number of iterations
#'
#' @importFrom is.circular circular
#'
#' @export
#' @useDynLib mvCircular
#'
#' @examples 
#' geomedian.circular(rmvvonmises(100,rep(0,3),rep(1,3),matrix(0,nrow=3,ncol=3)))
#' 
geomedian.circular <- function(data, tol = 1E-6, maxiter = 1E5) {
  
  
  if ( !is.mvCircular(data)) stop("All columns musth be circular")
  if ( !is.numeric(tol) || tol < 0 ) stop("Tolerance should be positive")
  if ( !is.numeric(maxiter) || maxiter < 0 || floor(maxiter) != maxiter ) stop("Maxiter should be a positive integer")
  
  # Get data size
  p <- ncol(data)
  n <- nrow(data)
  
  # Call procedure
  ret <- .C("_R__circularGeometricMedian",
            n = as.integer(n),
            p = as.integer(p),
            angles = as.double(t(data)),
            eps = as.double(tol),
            maxiter = as.integer(maxiter),
            maxiter_flag = integer(1),
            median = double(p),
            PACKAGE = "mvCircular")
  
  return( list(median = ret$median, maxiter_flag = ret$maxiter_flag) )
}