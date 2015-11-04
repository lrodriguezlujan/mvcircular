#' @name mvWrappedCauchy
#' @rdname mvWrappedCauchy
#' @title Multivariate wrapped Cauchy distribution
#' 
#' @description 
#' These functions implement diverse functionality over the 
#' multivariate wrapped normal distribution given its parameters mu, the circular mean vector,
#' and sigma, the variance covariance matrix.
#' 
#' @author Luis Rodriguez Lujan 
#' 
#' @keywords multivariate cauchy wrapped
#' 
#' @seealso \code{\link{mvCircularProbDist}}
#' @export
NULL

MVWRAPPEDCAUCHY_CLASS <- "mvWrappedCauchy"

#' @param mu Circular mean vector
#' @param sigma a positive-definite symmetric matrix that specifies covariance matrix
#' @param \dots (Constructor) Named list with additional attributes to add to the object
#' 
#' @examples 
#' mvWrappedCauchy(rep(0,3), diag(3) )
#' 
#' @importFrom circular is.circular conversion.circular as.circular
#' 
#' @export
mvWrappedCauchy <- function(mu, sigma, ...){
  
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
  
  # Add claseses (probdist + vonmises)
  class(obj) <- append( class(obj), MVWRAPPEDCAUCHY_CLASS)
  
  return(obj)
}


#' \code{mvWrappedCauchy.fit} uses data matrix or dataframe \code{samples} to compute the ML parameters of the distribution
#' 
#' @param samples Matrix or DF with mv circular samples
#' @param zero.threshold Any sigma value that verifies that \code{abs(x) < zero.threshold } is returned as zero
#' 
#' @return \code{mvWrappedCauchy.fit} returns a mvWrappedCacuhy object
#' 
#' @importFrom  Matrix nearPD
#' 
#' @examples 
#' samples <- rmvWrappedCauchy(1E5, rep(pi,3), matrix( c(0.3,.1,-.1,.1,.3,0,-.1,.0,.3), ncol = 3 , nrow = 3 )   )
#' obj <-mvWrappedCauchy.fit(samples)
#' sum(abs(obj$mu - rep(pi,3)))
#' sum(abs(obj$sigma - matrix( c(0.3,.1,-.1,.1,.3,0,-.1,.0,.3), ncol = 3 , nrow = 3 ) ))
#' plot(obj, data = obj$fitted.data[1:1000,])
#'
#' @rdname mvWrappedCauchy 
#' @export
mvWrappedCauchy.fit <- function(samples, zero.threshold = 1E-2){
  
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
  diag(sigma) <-  log( Mod( samples.complex.mean ) ) ^ 2
  
  # fill sigma off diagonal elements
  #
  for (j in 1:(ndim - 1) ) {
    for (k in (j + 1):ndim) {
      
      # Top
      aux <- Mod( mean( complex( argument = (samples[,j] + samples[, k]) ) ) )
      
      # Set matrix
      sigma[j, k] <- (log( aux ) ^ 2 - sigma[j, j] - sigma[k,k] ) / 2
      sigma[k, j] <- sigma[j, k]
      
    }
  }
  
  # Any value under the threshold is 0
  sigma[ abs(sigma) < zero.threshold ] <- 0
  
  # Sigma must be Positive Definite
  sigma <- as.matrix(Matrix::nearPD(sigma)$mat)
  
  # Return the object
  return( mvWrappedCauchy(mu, sigma, fitted.data = samples) )
}


#' @param n Number of samples to generate
#' @param ... (\code{rmvWrappedCauchy}) Additional parameters for \code{\link{mvtnorm::rmvnorm} }
#' 
#' @return \code{rmvWrappedCauchy} returns a multivariate circular dataframe with \code{n} 
#' samples from a wrapped-Cauchy distribution
#' 
#' @examples 
#' samples <- rmvWrappedCauchy(100, rep(pi,2), diag(2) )
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @importFrom corpcor is.positive.definite
#' @importFrom mvtnorm rmvnorm
#' 
#' @note 
#' \code{rmvWrappedCauchy} sampling  and density function method is based on the code from  the (discontinued) 
#' \code{LaplacesDaemon} package.
#' 
#' @rdname mvWrappedCauchy
#' @export
rmvWrappedCauchy <- function(n, mu, sigma, ...){
  
  # Check n
  if ( !is.numeric(n) || n <= 0 || floor(n) != n ) stop("The number of samples should be a positive integer")
  
  #FROM: LAPLACESDEMON
  
  # Mu as matrix
  mu <- rbind(mu)
  
  # Default values for sigma
  if (missing(sigma)) sigma <- diag( ncol(mu) )
  if (!is.matrix(sigma)) sigma <- matrix(sigma)
  
  # Check that sigma is PD
  if (!corpcor::is.positive.definite(sigma))
    stop("Matrix sigma is not positive-definite.")
  
  k <- ncol(sigma)
  if (n > nrow(mu)) mu <- matrix(mu, n, k, byrow = TRUE)
  
  # Sample from chisq
  x <- rchisq(n,1)
  x <- ifelse(x == 0, 1e-100, x)
  
  # Sample from multivariate normal
  z <- mvtnorm::rmvnorm(n, rep(0,k), sigma, ...)
  
  # Put in place
  x <- (mu + z/sqrt(x)) %% (2*pi)

  # Retun mv df
  return( as.mvCircular(x) )
}

# Auxiliar multivariate cauchy density (FROM: LAPLACESDEMON)
dmvc <- function(x, mu, sigma, log=FALSE)
{
  # Default values and checks
  if ( !is.matrix(x)) x <- rbind(x)
  if (missing(sigma)) sigma <- diag(ncol(x))
  if (!is.matrix(sigma)) sigma <- matrix(sigma)
  if (!corpcor::is.positive.definite(sigma))
    stop("Matrix S is not positive-definite.")
  
  k <- nrow(sigma)
  ss <- sweep(as.matrix(x),2, mu)
  omega <- MASS::ginv(sigma)
  z <- Matrix::rowSums({ss %*% omega} * ss)
  dens <- as.vector(lgamma(k/2) - (lgamma(0.5) + log(1 ^ (k/2)) +
                                     (k/2)*log(pi) + 0.5*log(det(sigma)) + ((1 + k)/2)*log(1 + z)))
  
  if (log == FALSE) dens <- exp(dens)
  
  return(dens)
}

#'  \code{dmvWrappedCauchy} computes multivariate wrapped normal densitiy function approximately. The precission is controled by
#'  \code{k}, the number of points to evaluate per dimension. The total number of points will be  \code{(k+1) ^ ndim} (Z-lattice)
#' 
#' @param x The point to evaluate
#' @param k Number of points per dimension
#' 
#' @examples
#' dmvWrappedCauchy(c(0,0,0), rep(0,3), 1000*diag(3) )
#' dmvWrappedCauchy(c(0,0,0), rep(0,3), 1000*diag(3),k = 20)
#' dmvWrappedCauchy(c(0,0,0), rep(0,3), 1000*diag(3),k = 2)
#' dmvWrappedCauchy(c(0,0,0), rep(0,3), 0.1*diag(3),k = 2)
#' dmvWrappedCauchy(c(pi,pi,pi), rep(0,3), diag(3) )
#' 
#' @importFrom corpcor is.positive.definite
#' @importFrom Matrix rowSums
#' @importFrom MASS ginv
#' 
#' @rdname mvWrappedCauchy
#' @export
dmvWrappedCauchy <- function(x, mu, sigma, k = 10){
  
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
    sum(dmvc( sweep(lattice, 2, y, FUN = "+" ),
              mu = mu, sigma = sigma, log  = F ))  
  },original.lat)
  
  # Return!
  return( as.numeric( val ) )
}

#' @examples 
#' obj <- mvWrappedCauchy(rep(pi,2), diag(2) )
#' samples <- getSamples(obj,100)
#' plot(as.numeric(samples[,1]), as.numeric(samples[,2]) )
#' 
#' @rdname mvWrappedCauchy
#' @export
getSamples.mvWrappedCauchy <- function(obj, n, ...) {
  
  # Retun mv df
  return( rmvWrappedCauchy(n, obj$mu, obj$sigma, ... ) )
}

#' @examples
#' obj <- mvWrappedCauchy(rep(0,3), diag(3) )
#' fval(obj,c(0,0,0))
#' fval(obj,c(2*pi,2*pi,2*pi))
#' fval(obj,c(pi,pi,pi))
#' 
#' obj <- mvWrappedCauchy(rep(0,3), 1000*diag(3) )
#' fval(obj,c(0,0,0), k= 60 ) # High accuracy
#' fval(obj,c(0,0,0), k= 20 )
#' fval(obj,c(0,0,0), k = 2 ) # Low accuracy due to dispersion
#' 
#' obj <- mvWrappedCauchy(rep(0,3), 0.1*diag(3) )
#' fval(obj,c(0,0,0), k= 60 ) # Low dispersion, almost all density is inside 0,2pi 
#' fval(obj,c(0,0,0), k= 2 ) # Low dispersion, almost all density is inside 0,2pi 
#' 
#' @rdname mvWrappedCauchy
#' @export
fval.mvWrappedCauchy <- function(obj, x, k = 10 ) {
  return(dmvWrappedCauchy(x, obj$mu, obj$sigma, k ))
}

#' @importFrom circular dwrappedcauchy
#' @importFrom MASS ginv
#' @rdname mvWrappedCauchy
circMarginal.mvWrappedCauchy <- function(obj, x, i){
  return( circular::dwrappedcauchy(x,obj$mu[i], exp(((obj$sigma)[i,i]) ^ 2 * -2 ) ) )
}

#' @rdname mvWrappedCauchy
circMarginalMean.mvWrappedCauchy <- function(obj , i){
  return( obj$mu[i] )
}

#' @importFrom MASS ginv
#' @rdname mvWrappedCauchy
circMarginalConcentration.mvWrappedCauchy <- function(obj, i){
  return( exp(((obj$sigma)[i,i]) ^ 2 * -2 ) )
}

#' @importFrom MASS ginv
#' @rdname mvWrappedCauchy
circCor.mvWrappedCauchy <- function(obj, i, j){
  return( MASS::ginv(obj$sigma)[i, j] )
}