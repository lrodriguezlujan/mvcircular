MVVONMISES_CLASS <- "mvVonMises"

#' Multivariate von Mises distribution
#' 
#' Creates a multivariate von Mises object
#' 
#' @importFrom circular is.circular conversion.circular
#' 
#' @param mu Circular mean vector, can be double
#' @param kappa Concentration vector
#' @param lambda 
#' 
#' @return a Multivariate von Mises S3 object
#' 
#' @examples 
#' mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' 
#' @export
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
  
  # Create list object
  obj <- list(
    dim = length(mu),
    mu = mu,
    kappa = kappa,
    lambda = lambda,
    z = NULL)
  
  ## Add additional params
  obj <- c(obj, list(...))
  
  # Add claseses (probdist + vonmises)
  class(obj) <- c(PROBDIST_CLASS, MVVONMISES_CLASS)
  
  return(obj)
}

#' Computes mv von Mises normalization term
#' 
#' Given a multivariate von Mises distribtuion, approximates the normalization
#' term using monte-carlo integration over a uniform mesh.
#' 
#' @param obj A multivariate von Mises distribution
#' @param normalization.samples Number of points to evaluate in the mc integration
#' 
#' @return Mv von mises object with  approximated normalization term set (z)
#' 
#' @examples 
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' obj <- normalize(obj, normalization.samples = 1E6)
#' obj$z
#' 
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

#' Fit von Mises distribution
#' 
#' This method uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with mv circular samples
#' @param ... Additional parameters for \see{fit_mvvonmises}
#' 
#' @return A mvVonMises object
#' 
#' @examples 
#' samples <- rmvVonMises(1E5, rep(0,4), rep(1,4), matrix(0,ncol=4,nrow=4) )
#' obj.fitted <- mvVonMises.fit(samples)
#' plot(obj, data = obj$fitted.data)
#' 
#' @export
mvVonMises.fit <- function(samples, ...){
  
  # Call fit_mvvonmises
  fit.ret <- fit_mvvonmises(samples, ...)
 
  # Check for error
  if (is.na(fit.ret$loss) )
    stop("Fitting procedure failed")
  
  # Create a mvVonMises
  return( mvVonMises(fit.ret$mu, fit.ret$kappa, fit.ret$lambda, 
                     fitted.data = samples, fitted.loss = fit.ret$loss, fitted.params = list(...) ) )
}

#' Multivariate von Mises sampler
#'
#' Samples n instances from the multivariate vonMises 
#'
#' @param obj Multivariate von Mises distribution
#' @param n Number of samples to generate
#' @param ... Additional parameters for \see{rmvvonmises_rs}
#' 
#' @return A circular df
#' 
#' @importFrom circular circular
#' 
#' @examples 
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' plot(getSamples(obj,100))
#' 
#' @export
getSamples.mvVonMises <- function(obj, n, ...) {
  
  # Check obj
  if ( !inherits(obj, MVVONMISES_CLASS) ) stop("Object should be a von Mises multivariate distribution")
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  # Call sampler
  res <- rmvvonmises_rs( n, obj$mu, obj$kappa, obj$lambda, ...)
  
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

#' ENGANCHAR A getSamples (ver doc S3)
#' @param n Number of samples to generate
#' @param mu Circular mean vector
#' @param kappa Concentration vector
#' @param lambda 
#' @param ... Additional parameters for \see{rmvvonmises_rs}
#' 
#' @export
rmvVonMises <- function(n, mu, kappa, lambda, ...){
  return( getSamples(mvVonMises(mu,kappa,lambda), n) )
}

#' Multivariate von Mises density function
#' 
#' Computes von Mises density function value. The normalization term can be computed via approximation (mc integration)
#' but result can be unaccurate
#' 
#' @param obj A Multivariate von Mises object
#' @param x The point to evaluate
#' 
#' @return Density function evaluated on the given point
#' 
#' @examples
#' obj <- mvVonMises(rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' # Compute normalization term
#' obj <- normalize(obj, normalization.samples = 1E6)
#' fval(obj, c(0,0,0))
#' dmvVonMises(c(0,0,0),rep(0,3), rep(1,3), matrix(0,ncol=3,nrow=3) )
#' 
#' @export
fval.mvVonMises <- function(obj, x, ... ) {
  
  # Validate inputs
  # Validate inputs
  if ( (is.matrix(x) || is.data.frame(x) ) ) {
    if (ncol(x) != length(obj$mu))  stop("")
  }
  else if (is.numeric(x)) {
    if ( length(x) != length(obj$mu) ) stop("")
    else x <- matrix(x,nrow = 1)
  }
  else stop("")
  
  # Compute unormalized value
  aux <- as.matrix(sweep(x,2,obj$mu)) # (x-mu) to each point
  
  # Get values
  val <- exp( as.numeric(tcrossprod(obj$kappa,cos(aux) )) 
              +  Matrix::rowSums(0.5 * ((sin(aux) %*% obj$lambda ) * sin(aux))) )
  
  if ( is.null(obj$z) ) {
    warning("Returning unnormalized value")
    return( as.numeric(val) )
  }
  else{
    return(  as.numeric(val)/obj$z)
  }
}

#'  ENGANCHAR A LA ANTERIOR
#' @param x 
#' @param mu Circular mean vector
#' @param kappa Concentration vector
#' @param lambda 
#' @param normalize
#' @param normalization.sampless
#' 
#' @export
dmvVonMises <- function(x, mu, kappa, lambda, normalization.samples = 1E6){
  obj <- mvVonMises(mu,kappa,lambda)
  if (normalization.samples > 0 ) obj <- normalize(obj,normalization.samples)
  return( fval(obj,x))
}

#' Multivariate von Mises plot
#' 
#' Plots a multivariate von mises distribution. The plot is a d x (d+1) grid with marginals on the diagonal, 
#' lambda coefficients in the upper tirangle and data plots on the lower.
#' 
#' @param obj Multivariate von mises distribution
#' @param data Datapoints to plot. Usually the samples used in \code{mvVonMises.fit}
#' @param n If data is null, number of point to sample from the distribution
#' @param ... Additional plot parameters \link{circularMvVm.plot}
#' 
#' @examples 
#' samples <- rmvVonMises(1E5, rep(0,4), rep(1,4), matrix(0,ncol=4,nrow=4) )
#' obj <- mvVonMises.fit(samples)
#' plot(obj, data = obj$fitted.data[1:100,])
#' 
#' 
# plot.mvVonMises <- function(obj, data = NULL, n = 1000, ...){
#   
#   # If data is null and n != 0 we create n samples 
#   if ( is.null(data) && n > 0 )
#     data <- getSamples(obj, n)
#   
#   # call plot_circularMvVm auxiliar func
#   circularMvVm.plot( obj$mu, obj$kappa, obj$lambda, data, ...)
# }

circMarginal.mvVonMises <- function(obj, x, i){
  return( circular::dvonmises(x,obj$mu[i], obj$kappa[i] ) )
}

circMarginalMean.mvVonMises <- function(obj , i){
  return( obj$mu[i] )
}

circMarginalConcentration.mvVonMises <- function(obj, i){
  return( obj$kappa[i] )
}

circCor.mvVonMises <- function(obj, i, j){
  return( obj$lambda[i, j] )
}