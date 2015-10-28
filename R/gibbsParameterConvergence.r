#' Vm Parameter inf. convergence
#'
#'
#'@export
gibbsParamVmConvergence <- function(n,mu,kappa,lambda,blocks=20,...){
  
  # Call method to Generate samples
  samples <- rmvvonmises_rs(n,mu,kappa,lambda,...)
  
  # Divide each one in blocks
  chainblocks <- lapply(samples,createBlocks,blocks)
  
  #Fit parameters for each chain and compute variance
  withinVariances <- Reduce(rbind,lapply(chainblocks,function(x){
    # X is a list of  blocks
    
    #Get parameters for each chain
    chainParams <- blockListParams(x)
    
    # Compute variances
    varianceVector <- vmParameterVariance(chainParams$mu,chainParams$kappa,chainParams$lambda)
    
    return(varianceVector)
  }))
  
  # Compute column average
  meanWithinVariance <- apply(withinVariances,2,mean)
  
  #
  # Mixed chain
  #
  
  # Compute mixed chain sample
  mixed <- mixMatricesAlternateRow(samples)
  mixedBlocks <- createBlocks(mixed,blocks)
  
  # Fit parameters for each block
  mixParams <- blockListParams(mixedBlocks)
  
  # Compute Variance vector
  mixVarVector <- vmParameterVariance(mixParams$mu,mixParams$kappa,mixParams$lambda)
  
  return(max(sqrt(meanWithinVariance/mixVarVector)))
}


#' Create data blocks
#'
#' Divide a sample in a list of n blocks
#'
#'
createBlocks <- function(data,nblocks){
  
  # Split in list
  p <- ncol(data)
  div <- split( t(data) , rep(1:nblocks, each = p*(nrow(data)/nblocks) ) )
  
  # tomatrix!
  return(lapply(div,function(x){matrix(x,ncol=p,byrow=T)}))
}

#' Mix dataset
#'
#' Mix data matrices alternating rows
#'
mixMatricesAlternateRow <- function(matrixlist){
  
  # Create final matrix
  p <- length(matrixlist)
  totalRows <- Reduce("+",lapply(matrixlist,nrow))
  res <- matrix(NA,ncol=ncol(matrixlist[[1]]),nrow=totalRows)
  
  # Copy to matrix
  for(i in 1:p){
    res[seq(i,totalRows,p),] <- matrixlist[[i]]
  }
  
  return(res)  
}

#' Vm Parameter varaiance
#' 
#' @param mu Matrix of mu parameters for each block
#' @param kappa Matrix of kappa parameters for each block
#' @param lambda Matrix of lambda parameters for each block
#'
#'@export
vmParameterVariance <- function(mu,kappa,lambda){
  
  # Compute mu circular variance
  mu_var<-apply(mu,2,var.circular)
  
  # For kappa and lambda apply linear variance
  kappa_var <- apply(kappa,2,stats::var)
  lambda_var <- apply(lambda,2,stats::var)
  
  # Return variance vector
  return(c(mu_var, kappa_var, lambda_var))
}

#' Compute block list parameters
#'
#'
#'
blockListParams <- function(blockList){
  
  #Get parameters for each block :) (in a list)
  chainParams <- lapply(blockList,function(y){
    params <- fit_mvvonmises(y)
    lambda <- matrix(params$lambda,ncol=ncol(y),nrow=ncol(y))
    return(list(mu=params$mu,kappa=params$kappa,lambda=lambda[lower.tri(lambda,diag = F)]))
  })
  
  # Put them together in matrices
  muchain <- Reduce(rbind,lapply(chainParams,function(x){return(x$mu)}))
  kappachain <- Reduce(rbind,lapply(chainParams,function(x){return(x$kappa)}))
  lambdachain <- Reduce(rbind,lapply(chainParams,function(x){return(x$lambda)}))
  
  return(list(mu=muchain,kappa=kappachain,lambda=lambdachain))
}