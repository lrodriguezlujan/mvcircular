#' Empiric KL divergence
#'
#' Calculates a KL divergence approximation using distribution samples
#' and a simplified circular distance
#'
#' @export
#' @useDynLib mvCircular
#'
#'@examples 
#' a <- rmvVonMises(10000,rep(0,3),rep(1,3),matrix(0,nrow=3,ncol=3))
#' b <- rmvVonMises(10000,rep(0,3),rep(1,3),matrix(0,nrow=3,ncol=3))
#' c <- rmvVonMises(10000,rep(2,3),rep(4,3),matrix(0,nrow=3,ncol=3))
#' div <- empKL.circular(a,b)
#' div$kl
#' div <- empKL.circular(a,c)
#' div$kl
empKL.circular <- function(a, b, n = nrow(a), m = nrow(b), k = 2, ...){
  
  # A can be a probability distribution
  if (inherits(a,"probDist")) {
    if (missing(n))
      stop("Needed number of samples to generate")
    else
      samples_a <- getSamples(a, n, ...)
  }
  else if ( is.mvCircular(a) ) 
    samples_a <- as.mvCircular(a)
  else 
    stop("")
  
  # Likewise for b
  if ( "probDist" %in% class(b) ) {
    if (missing(m))
      stop("Number of samples to generate should be given")
    else
      samples_b <- getSamples(b, n, ...)
  }
  else if ( is.mvCircular(a) ) 
    samples_b <- as.mvCircular(b)
  else
    stop("")
  
  # Now check that the dimensionality matches
  if (ncol(samples_a) != ncol(samples_b) )
    stop("Number of variables should match")
  
  # Next: Number of samples should be enough
  if (n < k || m < k) {
    warning("Not enough samples to compute empiric KL")
  }
  p <- ncol(samples_a)
  
  # FIXME
  if (k != 2) { 
    warning("Only implemented for k=2.")
    k <- 2
  }
  
  # Call procedure
  ret <- .C("__R_KL_estimator_circulardist",
            p = as.integer(p),
            k = as.integer(k),
            n = as.integer(n),
            A = as.double(t(samples_a)),
            m = as.integer(m),
            B = as.double(t(samples_b)),
            kl = double(1),
            PACKAGE = "mvCircular")

  # Get parameters
  return(list(sampleA = samples_a,
              sampleB = samples_b,
              k = k,
              kl = ret$kl))
}
