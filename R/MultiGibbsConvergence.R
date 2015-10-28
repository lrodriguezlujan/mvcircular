#' Multivariate gibbs convergence
#'
#' Computes convergence (PSRF) as defined by Brooks & Gelman 1998
#' 
#'
#'@export
MultiGibbsConvergence <- function(dataChains,meanFunction,distanceFunction,b=20){
  
  # TODO:
  return(0)
  
  # Aggregate all chains
  mixedChain <- reduce(rbind,dataChains)
  
  # Compute mixed mean
  mixedMean <- meanFunction(mixedChain)
  
  # Compute W (Intra variance)
  reduce("+",lapply(dataChains,function(x){
    # Compute mean
    chainMean <- meanFunction(x)
    # Compute component by component circular distance to the mean
    meanDistance <- apply(x,1,function(y){
      return(distanceFunction(chainMean,y)^2)
    })
  }))
  
  # Compute mean for every Chain
  chainMeans <- lapply(dataChain,meanFunction)

  # Distance to the mean for every chain
  chainDistances 
}