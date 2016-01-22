library(mvCircular)
library(ggplot2)
library(cowplot)
library(pracma)
library(prallel)

expset_ubias_A <- list()

### Function that executes the experiment
unbiasExp <- function(kappa, lambda, 
                   phi = NULL, #diag(kappa) - lambda,
                   H = NULL, #matrix(1, ncol = length(kappa), nrow = length(kappa)),
                   nrep = 1000, 
                   nmin = 10, 
                   nmax = 500, 
                   nfact = 2, 
                   cl = NULL){
  
  # Generate N values
  nlist <- nmin * nfact ^ (0:floor(log(nmax,nfact) - log(nmin,nfact)))
  # Add nmax if it is not in the list
  if (nlist[length(nlist)] != nmax) nlist <- c(nlist, nmax)
  
  ##
  # Create Distribution
  ##
  x <- mvVonMises( rep(0,length(kappa)), kappa = kappa, lambda)
  
  # Given parameter vector
  paramVector <- c(kappa, lambda[upper.tri(lambda)])
  
  # Execute experiment
  ret <- parallel::parLapplyLB(cl,nlist, function(n, nrep, real, phi, H, dist){
    #ret <- lapply(nlist, function(n, nrep, data, real){
    
    # Compute Penalized estimator
    vals <- sapply(1:nrep, function(x, dist, size){
      fit <- fit_mvvonmises( getSamples(dist, size, type = "rejection"), phi = phi, H = H)
      return(c(fit$kappa, fit$lambda[upper.tri(fit$lambda)]))
    },dist, n)
    
    # Transpose vals
    vals <- t(vals)
    
    # Expected Bias value
    bias <- base::colMeans(vals) - real
    
    # Return mean and variance
    return(list(n = n, bias = bias))
  },nrep = nrep, real = paramVector, phi = phi, H = H, dist = x )
  
  
  #results <- as.data.frame(t( sapply(ret, function(x){
  #  return( c( n = x$n, bias = x$bias) )
  #})))
  
  #colnames(results) <- c("n","bias")
  
  return(list(
    dist = x,
    results = ret
  #  results_table = results
  ))
}


################################
####### EXPERIMENT SETUP A #####
# 
# Eigenvalues range is small but with different centers (10,1,0.1)
# Lambda matrix is 0
#
################################

# CONFIG
nrep = 5000
nmax = 500
nmin = 10
nfact = 2
p = 3
ncores = 6

cl <- parallel::makeForkCluster(nnodes = ncores)

# Kappa with close-to-0 inf. norm
kappa_0 = rep(0.1,p)

# Kappa with greater than 1 inf. norm
kappa_10 = rep(1,p)

kappa_var = c(0.1,5,10)

lambda <- matrix(0, ncol = p, nrow = p)

## A.1
expset_ubias_A[[1]] <- unbiasExp(kappa_0, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

## A.2
expset_ubias_A[[2]] <- unbiasExp(kappa_10, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

expset_ubias_A[[3]] <- unbiasExp(kappa_var, lambda, NULL, NULL,
                                 nrep, nmin, nmax,
                                 nfact, cl)

parallel::stopCluster(cl)
