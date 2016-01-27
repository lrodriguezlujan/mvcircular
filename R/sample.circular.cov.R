#' Approximate circular covariance matrix
#'
#' @description Computes the sample circular covariance matrix
#'
#' @param d nxp data matrix
#' @param robust Use robust estimator for diagonal elements
#' 
#' @return pxp circular covariance matrix
#' 
#' @export
circular.cov <- function(d, robust = F) {
  
  if ( !is.mvCircular(d)) stop("Data should be multivariate circular")
  
  # Compute mean
  
  # First: Convert samples to complex numbers
  samples.complex <- data.frame( lapply(as.data.frame(d), function(x){ complex(argument = x) } ))
  
  # Then: compute mu component by component
  samples.complex.mean <- unlist(lapply(samples.complex, function(x){ mean(x) } ))
  
  # Mean vector
  mu <- Arg( samples.complex.mean ) %% (2*pi)
  
  # Rotate 
  d <- sweep(d,2, mu)
  d.s <- as.matrix(sin(d))
  cov.matrix <-  t(d.s) %*% d.s / nrow(d)
  
  if (!robust) {
    diag(cov.matrix) <- 2*(1 - colMeans(cos(d)))
  }
  
  return(cov.matrix)
}