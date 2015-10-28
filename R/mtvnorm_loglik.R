#' Computes negative log likelihood for a normal truncated distribution
#'
#'@export
#'@useDynLib mvCircular
#'
mvtnorm_loglik <- function(samples,mu = rep(0,ncol(samples)),
                           sigma = diag(ncol(samples)),
                           nSamples=1E5, lapack = T,
                           upper=rep(+Inf,ncol(samples)), lower=rep(-Inf,ncol(samples))){
  
  ## Approx norm term (montecarlo integration)
  mci_samples <- rmvnorm(nSamples,mean = mu,sigma = sigma)
  
  ## Samples within bounds
  F_val <- sum(apply(mci_samples,1,function(x,upp,low){
    return(x<=upp && x>=low)
  },upper,lower))/nSamples
  
  ## Centralize samples
  centralSamples <- sweep(samples,1,mu)
  
  if(!lapack){
    ## Sigma inverse
    s_inv <- ginv(sigma)
    partial<- centralSamples %*% s_inv
  }
  else{
    aux <- Matrix(sigma)
    c_sigma <- chol(aux)
    s_inv <- chol(c_sigma)
    partial <- centralSamples %*% s_inv
  }
  
  # Right part. sumprod 
  total <- 0
  for(i in 1:nrow(partial))
  {
    total <- total + sum(partial[i,] * centralSamples[i,])
  }
  
  return(nrow(samples)*(F_val + total)) 
}