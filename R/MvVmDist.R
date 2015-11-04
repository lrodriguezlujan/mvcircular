#' @name mvVonMises
#' @rdname mvVonMises
#' @title Multivariate von mises distribution
#' 
#' @description 
#' These functions implement diverse functionality over the 
#' multivariate von mises distribution given its parameters mu, the circular mean vector,
#' kappa, the concentration vector, and lambda, the dependency matrix.
#' 
#' @keywords multivariate vonmises
#' 
#' @author Luis Rodriguez Lujan 
#'  
#' @seealso \code{\link{mvCircularProbDist}}
#' @export
NULL

## Class name
MVVONMISES_CLASS <- "mvVonMises"

#' @importFrom circular is.circular conversion.circular as.circular
#' 
#' @param mu Circular mean vector, can be double
#' @param kappa Concentration vector
#' @param lambda Square depenedency matrix
#' @param \dots (\code{mvVonMises}) Named list with aditional parameters to include in the object
#' 
#' @examples 
#' mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' 
mvVonMises <- function(mu, kappa, lambda, ...){
  
  # Check parameters type
  if ( circular::is.circular(mu) ) mu <- circular::conversion.circular(mu,modulo = "2pi")
  else if ( is.numeric(mu) ) mu <- circular::as.circular(mu, modulo = "2pi", zero = 0, template = "none",
                                                         type = "angles", units = "radians", rotation = "counter" )
  
  else stop("Mu must be circular or numeric. see circular::as.circular() ")
  
  if ( !is.numeric(kappa) || any(kappa < 0) ) stop("Kappa should be a positive numeric vector")
  
  if ( !is.numeric(lambda) || !is.matrix(lambda) ) stop("Lambda should be a numeric matrix")
  
  # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
  if ( length(mu) != length(kappa) ) {stop("Parameters length do not match")}
  else if (length(mu) != nrow(lambda) ) {stop("Parameters length do not match")}
  else if (nrow(lambda) != ncol(lambda)) {stop("Lambda is not a square matrix")}
  
  # Create base object
  obj <- mvCircularProbDist( length(mu) , ...  )
  
  # Add params
  obj$mu <- mu
  obj$kappa <- kappa
  obj$lambda <- lambda

  # Add claseses (probdist + vonmises)
  class(obj) <- append( class(obj), MVVONMISES_CLASS)
  
  return(obj)
}



#' Fit method uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with mv circular samples
#' @param zero.threshold Any lambda value that verifies that \code{abs(x) < zero.threshold } is returned as zero
#' @param \dots (\code{fit}) Additional parameters for \code{\link{fit_mvvonmises}}
#'
#' 
#' @return \code{mvVonMises.fit} returns a multivariate VonMises object fitted from data
#' 
#' 
#' @examples 
#' samples <- rmvVonMises(1E5, rep(0,4), rep(1,4), matrix(0,ncol=4,nrow=4) )
#' obj.fitted <- mvVonMises.fit(samples)
#' plot(obj.fitted, data = obj.fitted$fitted.data[1:1000,])
#' 
#' @rdname mvVonMises
#' @export
mvVonMises.fit <- function(samples, zero.threshold = 1E-2, ...){
  
  # Call fit_mvvonmises
  fit.ret <- fit_mvvonmises(samples, ...)
 
  # Check for error
  if (is.na(fit.ret$loss) )
    stop("Fitting procedure failed")
  
  # Any value under the threshold is 0
  fit.ret$lambda[ abs(fit.ret$lambda) < zero.threshold ] <- 0
  
  # Create a mvVonMises
  return( mvVonMises(fit.ret$mu, fit.ret$kappa, fit.ret$lambda, 
                     fitted.data = samples, fitted.loss = fit.ret$loss, fitted.params = list(...) ) )
}



#' \code{rmvVonMises} samples n instances from a multivariate von Mises distribution
#' 
#' @param n Number of samples to generate
#' @param \dots (\code{rmvVonMises}) Additional parameters for \code{\link{rmvvonmises_rs}}
#' 
#' @return \code{rmvVonMises} returnrs a circular df with n samples
#' 
#' @importFrom circular circular
#' 
#' @rdname mvVonMises
#' @export
rmvVonMises <- function(n, mu, kappa, lambda, ...){
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  # Call sampler
  res <- rmvvonmises_rs( n, mu, kappa, lambda, ...)
  
  # rmvvonmises_rs samples from several chains (in gibbs mode). Join all values
  if (is.list(res)) {
    if ( length(res) == 1 )
      return( as.mvCircular( res[[1]] ) )
    else
      return( as.mvCircular( Reduce(rbind,res)[sample.int( n * length(res), size = n, replace = F), ] ))
  }
  else
    return( as.mvCircular(res  ) )
  

}

#' \code{dmvVonMises} computes von Mises density function value at the given points
#' 
#' @param x numeric or circular vector. Can also be a matrix
#' @param z normalization term
#' 
#' @importFrom Matrix rowSums
#' 
#' @rdname mvVonMises
#' @export
dmvVonMises <- function(x, mu, kappa, lambda, z = NULL){

  # Validate inputs
  if ( (is.matrix(x) || is.data.frame(x) ) ) {
    if (ncol(x) != length(mu))  stop("")
  }
  else if (is.numeric(x)) {
    if ( length(x) != length(mu) ) stop("")
    else x <- matrix(x,nrow = 1)
  }
  else stop("")
  
  # Compute unormalized value
  aux <- as.matrix(sweep(x,2,mu)) # (x-mu) to each point
  
  # Get values
  val <- exp( as.numeric(tcrossprod(kappa,cos(aux) )) 
              +  Matrix::rowSums(0.5 * ((sin(aux) %*% lambda ) * sin(aux))) )
  
  if ( is.null(z) ) {
    warning("Returning unnormalized value")
    return( as.numeric(val) )
  }
  else{
    return(  as.numeric(val)/z)
  }
  
}

#' \code{Normalize} approximates the normalization
#' term using monte-carlo integration over a uniform mesh.
#' 
#' @param obj A multivariate von Mises distribution
#' @param normalization.samples Number of points to evaluate in the mc integration
#' 
#' @return \code{normalize} returns a multivariate von Mises object with the approximated normalization term (z)
#' 
#' @examples 
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' obj <- normalize(obj, normalization.samples = 1E6)
#' obj$z
#' 
#' @rdname mvVonMises
#' @export
normalize.mvVonMises <- function(obj, normalization.samples = 1E6 ) {
  
  # Compute normalization term
  
  # Number of variables
  dims <- length(obj$mu)
  
  # Compute number point per dim
  samplesPerDim <- normalization.samples ^ (1/dims)
  
  # Integration hypervolume
  vol <- (2*pi) ^ dims
  
  # Generate samples matrix
  
  # Since 2pi == 0 we create  n+1 samples per dim and remove the last one (2pi) to avoid double-evaluating the same point
  # FIXME: We can take advantage of the symmetry of the distribution and "double" the performance
  dimvals <- seq(from = 0, to = 2*pi, length.out = (samplesPerDim + 1))[-(samplesPerDim + 1)] 
  points  <- expand.grid( split(rep(dimvals, each = dims),1:dims) )
  
  # Compute f value for each object 
  # NOTE: we can call .f function on each value, but that doenst seem really efficient.
  
  # Compute unormalized value
  aux <- as.matrix(sweep(points,2,obj$mu)) # (x-mu) to each point
  
  # Get values
  val <- exp( as.numeric(tcrossprod(obj$kappa,cos(aux) )) 
              +  Matrix::rowSums(0.5 * ((sin(aux) %*% obj$lambda ) * sin(aux))) )
  
  # Compute value
  obj$z <- vol/normalization.samples * sum(val) ;
  
  return(obj)
}

#' @examples 
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' plot(getSamples(obj,100))
#' 
#' @rdname mvVonMises
#' @export
getSamples.mvVonMises <- function(obj, n, ...) {
  
  # Check obj
  if ( !inherits(obj, MVVONMISES_CLASS) ) stop("Object should be a von Mises multivariate distribution")
  return( rmvVonMises( n, obj$mu, obj$kappa, obj$lambda ) )
}

#' @examples
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' # Compute normalization term
#' obj <- normalize(obj, normalization.samples = 1E6)
#' fval(obj, c(0,0,0))
#' # unnormalzied value
#' dmvVonMises(c(0,0,0),rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' 
#' @rdname mvVonMises
#' @export
fval.mvVonMises <- function(obj, x, ... ) {
  if ( !inherits(obj, MVVONMISES_CLASS) ) stop("Object should be a von Mises multivariate distribution")
  return( dmvVonMises(x, obj$mu, obj$kappa, obj$lambda, obj$z) )
}

#' @importFrom circular dvonmises
#' @rdname mvVonMises
#' @export
circMarginal.mvVonMises <- function(obj, x, i){
  return( circular::dvonmises(x,obj$mu[i], obj$kappa[i] ) )
}

#' @rdname mvVonMises
#' @export
circMarginalMean.mvVonMises <- function(obj , i){
  return( obj$mu[i] )
}

#' @rdname mvVonMises
#' @export
circMarginalConcentration.mvVonMises <- function(obj, i){
  return( obj$kappa[i] )
}

#' @rdname mvVonMises
#' @export
circCor.mvVonMises <- function(obj, i, j){
  return( obj$lambda[i, j] )
}