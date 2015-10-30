MVWRAPPEDNORMAL_CLASS <- "mvWrappedNormal"

#' Multivariate Wrapped normal distribution
#' 
#' Creates a multivariate wrapped normal object
#' 
#' 
#' @param mu Circular mean vector
#' @param sigma a positive-definite symmetric matrix that specifies covariance matrix
#' 
#' @return a Multivariate wrapped normal S3 object
#' 
#' @examples 
#' mvWrappedNormal(rep(0,3), diag(3) )
#' 
#' @export
mvWrappedNormal <- function(mu, sigma, ...){
  
  # Check parameters type
  if ( circular::is.circular(mu) ) mu <- circular::conversion.circular(mu,modulo = "2pi")
  else if ( is.numeric(mu) ) mu <- circular::as.circular(mu, modulo = "2pi", zero = 0, template = "none",
                                                         type = "angles", units = "radians", rotation = "counter" )
  
  else stop("Mu must be circular or numeric. see circular::as.circular() ")
  
  if ( !is.numeric(sigma) || !is.matrix(sigma) ) stop("Lambda should be a numeric matrix")
  
  # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
  else if (length(mu) != nrow(sigma) ) {stop("Parameters length do not match")}
  else if (nrow(sigma) != ncol(sigma)) {stop("sigma is not a square matrix")}
  
  # Create list object
  obj <- list(
    mu = mu,
    sigma = sigma)
  
  ## Add additional params
  obj <- c(obj, list(...))
  
  # Add claseses (probdist + vonmises)
  class(obj) <- c(PROBDIST_CLASS, MVWRAPPEDNORMAL_CLASS)
  
  return(obj)
}


#' Fit wrapped normal distribution
#' 
#' This method uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with mv circular samples
#' 
#' @return A mvWrappedNormal object
#' 
#' @examples 
#' samples <- rmvWrappedNormal(100000, rep(pi,3), matrix( c(3,1,-1,1,3,0,-1,0,3), ncol = 3 , nrow = 3 )   )
#' obj <- mvWrappedNormal.fit(samples)
#' sum(abs(obj$mu - rep(pi,3)))
#' sum(abs(obj$sigma - matrix( c(3,1,-1,1,3,0,-1,0,3), ncol = 3 , nrow = 3 ) ))
#' plot(obj)
#' 
#' @export
mvWrappedNormal.fit <- function(samples, zero.threshold = 1E-2, ...){
  
  # number of variables
  ndim <- ncol(samples)
  nsamples <- nrow(samples)
  
  # Prealloc sigma.sinh
  sigma <- matrix(numeric(1), ncol = ndim, nrow = ndim )
  
  # First: Convert samples to complex numbers
  samples.complex <- data.frame( lapply(samples, function(x){ complex(argument = x) } ))
  
  # Then: compute mu component by component
  samples.complex.mean <- unlist(lapply(samples.complex, function(x){ mean(x) } ))
  
  # Mean vector
  mu <- Arg( samples.complex.mean )
  
  # Sigma diagonal elements
  diag(sigma) <- -2 * log( Mod( samples.complex.mean ) )
  
  
  # fill sigma off diagonal elements
  #
  for (j in 1:(ndim - 1) ) {
    for (k in (j + 1):ndim) {
      
      # Top
      top <- Mod( mean( complex( argument = (samples[,j] - samples[, k]) ) ) )
      
      # Bottom
      bottom <- Mod( mean( complex( argument = (samples[,j] + samples[, k]) ) ) )
      
      # Set matrix
      sigma[j, k] <- log(top / bottom) / 2
      sigma[k, j] <- sigma[j, k]
    }
  }
  
  # Any value under the threshold is 0
  sigma[ abs(sigma) < zero.threshold ] <- 0
  
  # Sigma must be Positive Definite
  sigma <- as.matrix(Matrix::nearPD(sigma)$mat)
  
  # Return the object
  return( mvWrappedNormal(mu, sigma, fitted.data = samples) )
}

#' Multivariate wrapped normals sampler
#'
#' Samples n instances from the multivariate wrapped normal
#'
#' @param obj Multivariate wrapped normal
#' @param n Number of samples to generate
#' @param ... Additional parameters for \link{mvtnorm::rmvnorm}
#' 
#' @return A circular dataframe
#' 
#' @importFrom mvtnorm rmvnorm
#' 
#' @examples 
#' obj <- mvWrappedNormal(rep(pi,2), diag(2) )
#' samples <- getSamples(obj,100)
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @export
getSamples.mvWrappedNormal <- function(obj, n, ...) {
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  # Call sampler
  res <- mvtnorm::rmvnorm(n, mean = obj$mu, sigma = obj$sigma, ...) %% (2*pi)
  
  # Retun mv df
  return( as.mvCircular(res  ) )
}

#' ENGANCHAR A getSamples (ver doc S3)
#' @param n Number of samples to generate
#' @param mu Circular mean vector
#' @param sigma
#' @param ... Additional parameters for \see{mvtnorm::rmvnorm}
#' 
#' @examples 
#' samples <- rmvWrappedNormal(100, rep(pi,2), diag(2) )
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @export
rmvWrappedNormal <- function(n, mu, sigma, ...){
  return( getSamples(mvWrappedNormal(mu, sigma), n, ...) )
}

#' Multivariate wrapped normal density function
#' 
#' Computes multivariate wrapped normal densitiy function approximately. The precission is controled by
#' k, the number of points to evaluate per dimension. The total number of points will be  (k+1) ^ ndim (Z-lattice)
#' 
#' @param obj A Multivariate wrapped normal object
#' @param x The point to evaluate
#' @param k Number of terms to be used per dimension
#' @param ... Extra arguments for \link{mvtnorm::dmvnorm}
#' 
#' @return Density function evaluated on the given point
#' 
#' @importFrom mvtnorm dmvnorm
#' 
#' @examples
#' obj <- mvWrappedNormal(rep(0,3), diag(3) )
#' fval(obj,c(0,0,0))
#' fval(obj,c(2*pi,2*pi,2*pi))
#' fval(obj,c(pi,pi,pi))
#' 
#' obj <- mvWrappedNormal(rep(0,3), 1000*diag(3) )
#' fval(obj,c(0,0,0), k= 60 ) # High accuracy
#' fval(obj,c(0,0,0), k= 20 )
#' fval(obj,c(0,0,0), k = 2 ) # Low accuracy due to dispersion
#' 
#' obj <- mvWrappedNormal(rep(0,3), 0.1*diag(3) )
#' fval(obj,c(0,0,0), k= 60 ) # Low dispersion, almost all density is inside 0,2pi 
#' fval(obj,c(0,0,0), k= 2 ) # Low dispersion, almost all density is inside 0,2pi 
#' @export
fval.mvWrappedNormal <- function(obj, x, k = 10 ) {
  
  # Validate inputs
  if ( (is.matrix(x) || is.data.frame(x) ) ) {
    if (ncol(x) != length(obj$mu))  stop("")
  }
  else if (is.numeric(x)) {
    if ( length(x) != length(obj$mu) ) stop("")
    else x <- matrix(x,nrow = 1)
  }
  else stop("")
  
  # Number of variables
  dims <- ncol(x)
  
  # Put x in 0,2pi, so computation is more accurate (place it where most prob is)
  x <- x %% (2*pi)
  
  # Make k even
  if ( (k %% 2) != 0) k <- k + 1
  
  
  # Symmetric values around 0
  dimvals <- seq( from = 2*pi, by = 2*pi, length.out = (k/2)  )
  dimvals <- c(-dimvals, 0, dimvals)
  
  # Create lattice
  original.lat <- expand.grid( split(rep(dimvals, each = dims),1:dims) )
  
  # Compute one time per point. We could have created a huge lattice, but here memory/cpu tradeoff isnt worth it
  val <- apply(x,1, function(y,lattice){
    # Evaluate at the given points
    # Points to evaluate (x +- k2pi)
    sum(mvtnorm::dmvnorm( sweep(lattice, 2, y, FUN = "+" ),
                                mean = obj$mu, sigma = obj$sigma, log  = F ))  
  },original.lat)
  
  # Return!
  return( as.numeric( val ) )
}

#'  ENGANCHAR A LA ANTERIOR
#' @param x 
#' @param mu Circular mean vector
#' @param sigma
#' @param k
#' @param ...
#' 
#' @export
#' 
#' @examples
#' dmvWrappedNormal(c(0,0,0), rep(0,3), 1000*diag(3) )
#' dmvWrappedNormal(c(0,0,0), rep(0,3), 1000*diag(3),k = 20)
#' dmvWrappedNormal(c(0,0,0), rep(0,3), 1000*diag(3),k = 2)
#' dmvWrappedNormal(c(0,0,0), rep(0,3), 0.1*diag(3),k = 2)
#' dmvWrappedNormal(c(pi,pi,pi), rep(0,3), diag(3) )
dmvWrappedNormal <- function(x, mu, sigma, k = 10, ...){
  obj <- mvWrappedNormal(mu,sigma)
  return( fval(obj, x, k, ...))
}

#' Multivariate wrapped normal plot
#' 
#' Plots a multivariate wrapped normal distribution. The plot is a d x (d+1) grid with marginals on the diagonal, 
#' sigma coefficients in the upper tirangle and data plots on the lower.
#' 
#' @param obj Multivariate wrapped normal distribution
#' @param data Datapoints to plot. Usually the samples from \code{mvWrappedNormal.fit}
#' @param n If data is null, number of point to sample from the distribution
#' @param ... Additional plot parameters \link{circularNorm.plot}
#' 
#' @examples 
#' samples <- rmvWrappedNormal(1000, rep(pi,3), matrix( c(0.3,0.1,-0.1,0.1,0.3,0,-0.1,0,0.3), ncol = 3 , nrow = 3 )   )
#' obj <- mvWrappedNormal.fit(samples)
#' plot(obj, data = obj$fitted.data )
#' 
#' 
#' @export
plot.mvWrappedNormal <- function(obj, data = NULL, n = 1000, ...){
  
  # If data is null and n != 0 we create n samples 
  if ( is.null(data) && n > 0 )
    data <- getSamples(obj, n)
  
  # call plot_circularMvVm auxiliar func
  wrappedDist.plot( obj$mu, obj$sigma, data, circular::dwrappednormal)
}