#' @name mvCircTruncNormal
#' @rdname mvCircTruncNormal
#' @title Multivariate circular truncated normal distribution
#' 
#' @description 
#' These functions implement diverse functionality over the 
#' multivariate wrapped normal distribution given its parameters mu, the circular mean vector,
#' and sigma, the variance covariance matrix.
#' 
#' @author Luis Rodriguez Lujan 
#' 
#' @keywords multivariate normal circular truncated
#' 
#' @seealso \code{\link{mvCircTruncNormal}}
#' @export
NULL

MVCIRCTRUNCNORMAL_CLASS <- "mvCircTruncNormal"


#' @param mu Circular mean vector
#' @param sigma a positive-definite symmetric matrix that specifies covariance matrix
#' @param \dots (Constructor) Named list with additional attributes to add to the object
#' 
#' @examples 
#' mvCircTruncNormal(rep(0,3), diag(3) )
#' 
#' @importFrom circular is.circular conversion.circular as.circular
#' 
#' @export
mvCircTruncNormal <- function(mu, sigma, ...){
  
  # Check parameters type
  if ( circular::is.circular(mu) ) mu <- circular::conversion.circular(mu,modulo = "2pi")
  else if ( is.numeric(mu) ) mu <- circular::as.circular(mu, modulo = "2pi", zero = 0, template = "none",
                                                         type = "angles", units = "radians", rotation = "counter" )
  
  else stop("Mu must be circular or numeric. see circular::as.circular() ")
  
  if ( !is.numeric(sigma) || !is.matrix(sigma) ) stop("Lambda should be a numeric matrix")
  
  # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
  else if (length(mu) != nrow(sigma) ) {stop("Parameters length do not match")}
  else if (nrow(sigma) != ncol(sigma)) {stop("sigma is not a square matrix")}
  
  # Create base object
  obj <- mvCircularProbDist( length(mu), ... )
  
  # Create list object
  obj$mu <- mu
  obj$sigma <- sigma
  obj$lower = as.numeric(mu) - pi
  obj$upper = as.numeric(mu) + pi
  
  # Add claseses (probdist + vonmises)
  class(obj) <- append( class(obj), MVCIRCTRUNCNORMAL_CLASS)
  
  return(obj)
}


#' Fit method uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with multivariate circular samples
#' @param zero.threshold Any sigma value that verifies that \code{abs(x) < zero.threshold } is returned as zero
#' 
#' @return \code{mvCircTruncNormal} returns a mvCircTruncNormal object
#' 
#' @importFrom  Matrix nearPD
#' @importFrom tmvtnorm mle.tmvnorm
#' 
#' @examples 
#' samples <- rmvCircTruncNormal(100000, rep(pi,3), matrix( c(3,1,-1,1,3,0,-1,0,3), ncol = 3 , nrow = 3 )   )
#' obj <- mvCircTruncNormal.fit(samples)
#' sum(abs(obj$mu - rep(pi,3)))
#' sum(abs(obj$sigma - matrix( c(3,1,-1,1,3,0,-1,0,3), ncol = 3 , nrow = 3 ) ))
#' plot(obj, obj$fitted.data [1:1000,])
#' 
#' @rdname mvCircTruncNormal
#' @export
mvCircTruncNormal.fit <- function(samples, zero.threshold = 1E-2, ...){
  
  # number of variables
  ndim <- ncol(samples)
  nsamples <- nrow(samples)
  
  # First: Convert samples to complex numbers
  samples.complex <- data.frame( lapply(samples, function(x){ complex(argument = x) } ))
  
  # Then: compute mu component by component
  samples.complex.mean <- unlist(lapply(samples.complex, function(x){ mean(x) } ))
  
  # Mean vector lowr and upper bounds
  mu <- Arg( samples.complex.mean ) %% (2*pi)
  lower <- as.numeric(mu) - pi
  upper <- as.numeric(mu) + pi
  
  # Put points in the region
  
  samples <- sweep( (samples %% (2*pi)) - pi, 2, mu, FUN = "+" )
  
  fit <- tmvtnorm::mle.tmvnorm(
    as.matrix( samples ),
    lower = lower,
    upper = upper,
    start = list( mu = as.numeric(mu), sigma = diag(ncol(samples) ) ),
    fixed = list( "mu" ))
  
  # Sigma elements
  sigma <- matrix(0, ncol = ndim, nrow = ndim)
  sigma[lower.tri(sigma,diag = T)] <- attr(fit,"coef")[-(1:ndim)]
  sigma <- sigma + t(sigma)
  diag(sigma) <- diag(sigma)/2
  
  # Any value under the threshold is 0
  sigma[ abs(sigma) < zero.threshold ] <- 0
  
  # Sigma must be Positive Definite (dk if mle method does it)
  sigma <- as.matrix(Matrix::nearPD(sigma)$mat)
  
  # Return the object
  return( mvCircTruncNormal(mu, sigma, fitted.data = samples) )
}




#' @param n Number of samples to generate
#' @param ... (\code{rmvCircTruncNormal}) Additional parameters for \code{\link{tmvtnorm::rtmvnorm} }
#' 
#' @return \code{rmvCircTruncNormal} returns a multivariate circular dataframe with \code{n} 
#' samples from a truncated circular normal distribution
#' 
#' @examples 
#' samples <- rmvCircTruncNormal(100, rep(pi,2), diag(2) )
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @importFrom tmvtnorm rtmvnorm
#' 
#' @rdname mvCircTruncNormal
#' @export
rmvCircTruncNormal <- function(n, mu, sigma, ...){
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  mu <- as.numeric(mu) %% (2*pi)
  
  # Call sampler
  res <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = sigma, lower =  mu - pi, upper = mu + pi, ...) %% (2*pi)
  
  # Retun mv df
  return( as.mvCircular(res  ) )
}



#' \code{dmvCircTruncNormal} computes multivariate circular normal densitiy function approximately. 
#' 
#' @param x The point to evaluate
#' @param ...  (\code{dmvCircTruncNormal}) extra arguments for \code{\link{tmvtnorm::dtmvnorm}}
#' 
#' @return \code{dmvCircTruncNormal} returns the density function evaluated at \code{x}
#' 
#' @importFrom tmvtnorm dtmvnorm
#' 
#' @rdname mvCircTruncNormal
#' @export
#' 
#' @examples
#' dmvCircTruncNormal(c(0,0,0), rep(0,3), 1000*diag(3) )
#' dmvCircTruncNormal(c(0,0,0), rep(0,3), 0.1*diag(3))
#' dmvCircTruncNormal(c(pi,pi,pi), rep(0,3), diag(3) )
dmvCircTruncNormal <- function(x, mu, sigma, ...){
  
  # Validate inputs
  if ( (is.matrix(x) || is.data.frame(x) ) ) {
    if (ncol(x) != length(mu))  stop("")
  }
  else if (is.numeric(x)) {
    if ( length(x) != length(mu) ) stop("")
    else x <- matrix(x,nrow = 1)
  }
  else stop("")
  
  # Number of variables
  dims <- ncol(x)
  
  
  # Put x in mu-pi,mu+pi
  mu <- as.numeric(mu) %% (2*pi)
  x <- (x %% (2*pi) ) + mu - pi
  
  return( tmvtnorm::dtmvnorm(x, 
                     mean = mu,
                     sigma = sigma,
                     lower = mu - pi,
                     upper = mu + pi) )
}

#' @examples 
#' obj <- mvCircTruncNormal(rep(pi,2), diag(2) )
#' samples <- getSamples(obj,100)
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @rdname mvCircTruncNormal
#' @export
getSamples.mvCircTruncNormal <- function(obj, n, ...) {
  # Retun mv df
  return( rmvCircTruncNormal(n, obj$mu, obj$sigma  ) )
}


#' @examples
#' obj <- mvCircTruncNormal(rep(0,3), diag(3) )
#' fval(obj,c(0,0,0))
#' fval(obj,c(2*pi,2*pi,2*pi))
#' fval(obj,c(pi,pi,pi))
#' 
#' obj <- mvCircTruncNormal(rep(0,3), 1000*diag(3) )
#' fval(obj,c(0,0,0))
#' 
#' obj <- mvCircTruncNormal(rep(0,3), 0.1*diag(3) )
#' fval(obj,c(0,0,0))
#' 
#' @rdname mvCircTruncNormal
#' @export
fval.mvCircTruncNormal <- function(obj, x, ... ) {
  return( dmvCircTruncNormal(x, obj$mu, obj$sigma, ...) )
}

#'@importFrom tmvtnorm dtmvnorm.marginal
#'@rdname mvCircTruncNormal
circMarginal.mvCircTruncNormal <- function(obj, x, i){
  return( tmvtnorm::dtmvnorm.marginal(x, n = i, mean = as.numeric(obj$mu),
                                      sigma = obj$sigma, lower = obj$lower, upper = obj$upper ) )
}

#'@rdname mvCircTruncNormal
circMarginalMean.mvCircTruncNormal <- function(obj , i){
  return( obj$mu[i] )
}

#'@rdname mvCircTruncNormal
circMarginalConcentration.mvCircTruncNormal <- function(obj, i){
  return( exp(((obj$sigma)[i,i]) ^ 2 / -2 ) )
}

#'@importFrom MASS ginv
#'@rdname mvCircTruncNormal
circCor.mvCircTruncNormal <- function(obj, i, j){
  return( MASS::ginv(obj$sigma)[i,j] )
}