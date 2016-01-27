#' @name mvNormal
#' @rdname mvNormal
#' @title Multivariate normal distribution

#' @author Luis Rodriguez Lujan 
#' 
#' @keywords multivariate normal 
#' 
#' @seealso \code{\link{mvCircularProbDist}}
#' @export
NULL

MVNORMAL_CLASS <- "mvNormal"


#' @param mu Circular mean vector
#' @param sigma a positive-definite symmetric matrix that specifies covariance matrix
#' @param \dots (Constructor) Named list with additional attributes to add to the object
#' 
#' @examples 
#' mvNormal(rep(0,3), diag(3) )
#' 
#' 
#' @export
mvNormal <- function(mu, sigma, ...){
  
  # Check parameters type
  if ( !is.numeric(mu) ) mu <- as.numeric(mu)
  
  if ( !is.numeric(sigma) || !is.matrix(sigma) ) stop("Lambda should be a numeric matrix")
  
  # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
  else if (length(mu) != nrow(sigma) ) {stop("Parameters length do not match")}
  else if (nrow(sigma) != ncol(sigma)) {stop("sigma is not a square matrix")}
  
  # Create base object
  obj <- mvCircularProbDist( length(mu), ... )
  
  # Create list object
  obj$mu <- mu
  obj$sigma <- sigma
  
  # Add claseses (probdist + vonmises)
  class(obj) <- append( class(obj), MVNORMAL_CLASS)
  
  return(obj)
}


#' Fit method uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with multivariate samples
#' @param zero.threshold Any sigma value that verifies that \code{abs(x) < zero.threshold } is returned as zero
#' 
#' @return \code{mvNormal} returns a mvNormal object
#' 
#' @importFrom  Matrix nearPD
#' 
#' 
#' @rdname mvNormal
#' @export
mvNormal.fit <- function(samples, zero.threshold = 1E-2, ...){
  
  # number of variables
  ndim <- ncol(samples)
  nsamples <- nrow(samples)
  
  # Prealloc sigma.sinh
  sigma <- stats::var(samples)
  mu <- colMeans(samples)
  
  # Any value under the threshold is 0
  sigma[ abs(sigma) < zero.threshold ] <- 0
  
  # Sigma must be Positive Definite
  sigma <- as.matrix(Matrix::nearPD(sigma)$mat)
  
  # Return the object
  return( mvNormal(mu, sigma, fitted.data = samples) )
}

#' \code{dmvNormal} computes multivariate normal densitiy function approximately. The precission is controled by
#' \code{k}, the number of points to evaluate per dimension. The total number of points will be  (k+1) ^ ndim (Z-lattice)
#' 
#' @param x The point to evaluate
#' @param ...  (\code{dmvNormal}) extra arguments for \code{\link{mvtnorm::dmvnorm}}
#' 
#' @return \code{dmvNormal} returns the density function evaluated at \code{x}
#' 
#' @importFrom mvtnorm dmvnorm
#' 
#' @rdname mvNormal
#' @export

dmvNormal <- function(x, mu, sigma, ...){
  
  # Validate inputs
  if ( (is.matrix(x) || is.data.frame(x) ) ) {
    if (ncol(x) != length(mu))  stop("")
  }
  else if (is.numeric(x)) {
    if ( length(x) != length(mu) ) stop("")
    else x <- matrix(x,nrow = 1)
  }
  else stop("")
  
  return(mvtnorm::dmvnorm( x, mean = mu, sigma = sigma, log  = F ))
  
}

#' @param n Number of samples to generate
#' @param ... (\code{rmvNormal}) Additional parameters for \code{\link{mvtnorm::rmvnorm} }
#' 
#' @return \code{rmvNormal} returns a multivariate circular dataframe with \code{n} 
#' samples from a normal distribution
#' 
#' @rdname mvNormal
#' @export
rmvNormal <- function(n, mu, sigma, ...){
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  # Call sampler
  res <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma, ...)
  
  # Retun mv df
  return( res )
}

#' 
#' @rdname mvNormal
#' @export
getSamples.mvNormal <- function(obj, n, ...) {
  # Retun mv df
  return( rmvNormal(n, obj$mu, obj$sigma  ) )
}


#' @rdname mvNormal
#' @export
fval.mvNormal <- function(obj, x, k = 10, ... ) {
  return( dmvNormal(x, obj$mu, obj$sigma, k , ...) )
}

#'@rdname mvNormal
#'@export
circMarginal.mvNormal <- function(obj, x, i){
  return( stats::dnorm(x,obj$mu[i], sd = obj$sigma[i,i] ) )
}

#'@rdname mvNormal
#'#'@export
circMarginalMean.mvNormal <- function(obj , i){
  return( obj$mu[i] )
}